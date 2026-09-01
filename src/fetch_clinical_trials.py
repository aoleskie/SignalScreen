#!/usr/bin/env python3
import requests
import pandas as pd
from pathlib import Path
import time
import re

from src.compound_names import is_searchable_term, strip_source_tag

CT_BASE = "https://clinicaltrials.gov/api/v2/studies"


def _normalize(text):
    """Lowercase, collapse punctuation to spaces, for tolerant name matching."""
    return re.sub(r"[^a-z0-9]+", " ", (text or "").lower()).strip()


def _tokens(text):
    return _normalize(text).split()


# Words that may legitimately precede a drug name in an intervention label:
# route, formulation and dosing qualifiers. Anything else before the name means
# it is a different chemical entity -- '13-cis-retinoic acid' is not
# 'retinoic acid', and 'isotretinoin' is not 'tretinoin'.
_ALLOWED_PREFIX_WORDS = {
    "oral", "topical", "intravenous", "iv", "injection", "injectable",
    "inhaled", "subcutaneous", "systemic", "low", "high", "dose", "doses",
    "single", "daily", "twice", "once", "micronized", "liposomal", "liposome",
    "cream", "gel", "lotion", "ointment", "foam", "solution", "capsule",
    "capsules", "tablet", "tablets", "sustained", "release", "extended",
    "placebo", "matching", "drug", "study", "arm", "group", "the", "a", "of",
}


def intervention_names(study: dict):
    """Every intervention name and alias declared by a study."""
    module = (study.get("protocolSection", {}) or {}).get("armsInterventionsModule", {}) or {}
    names = []
    for iv in module.get("interventions", []) or []:
        if iv.get("name"):
            names.append(iv["name"])
        names.extend(iv.get("otherNames") or [])
    return names


def study_matches_term(study: dict, term: str):
    """
    True when the compound is actually an intervention in this study.

    ClinicalTrials.gov matches loosely: query.intr='Vitamin A acid' returns 1888
    studies, nearly all of which are about some other vitamin. Requiring the term
    to appear in a declared intervention name is what makes the count mean
    something. Returns the matching intervention name for provenance.
    """
    needle = _tokens(term)
    if not needle:
        return None

    for name in intervention_names(study):
        hay = _tokens(name)
        for i in range(len(hay) - len(needle) + 1):
            if hay[i:i + len(needle)] != needle:
                continue
            # Anchored at the start, or preceded only by route/formulation words.
            if i == 0 or all(w in _ALLOWED_PREFIX_WORDS for w in hay[:i]):
                return name
    return None


# ------------------------------------------------------------
# Low-level API call
# ------------------------------------------------------------
def query_clinical_trials(term: str, max_studies: int = 200) -> list[dict]:
    """
    Call ClinicalTrials.gov v2 API using query.term=<term>
    Returns the list of study objects.
    """
    if not term:
        return []

    # query.intr searches intervention fields specifically. query.term searches
    # everything, which is how a compound picked up trials that merely mention a
    # related word somewhere in the protocol.
    params = {
        "query.intr": term,
        "pageSize": max_studies,
    }
    r = requests.get(CT_BASE, params=params, timeout=30)
    try:
        r.raise_for_status()
        data = r.json()
        return data.get("studies", [])
    except Exception as e:
        print(f"[query_clinical_trials] Error for term '{term}': {e}")
        return []


# ------------------------------------------------------------
# Correct parser (based on your example.json structure)
# ------------------------------------------------------------
def parse_clinical_study(study: dict) -> dict:
    """
    Extract required fields from a ClinicalTrials.gov v2 'study' object.

    The correct structure (from your example):
      - study["protocolSection"]["designModule"]["phases"]
      - study["protocolSection"]["conditionsModule"]["conditions"]
      - study["protocolSection"]["sponsorCollaboratorsModule"]["leadSponsor"]["name"]
      - study["protocolSection"]["statusModule"]["overallStatus"]
      - study["protocolSection"]["identificationModule"]["nctId"]
      - study["hasResults"]
    """

    protocol = study.get("protocolSection", {}) or {}

    # --- NCT ID ---
    ident = protocol.get("identificationModule", {}) or {}
    nct_id = ident.get("nctId")

    # --- Status ---
    status_mod = protocol.get("statusModule", {}) or {}
    status = status_mod.get("overallStatus")

    # --- Phase ---
    design_mod = protocol.get("designModule", {}) or {}
    phases = design_mod.get("phases") or []
    if isinstance(phases, list) and phases:
        phase = "; ".join(phases)
    else:
        phase = None

    # --- Conditions / Indication ---
    cond_mod = protocol.get("conditionsModule", {}) or {}
    conditions = cond_mod.get("conditions") or []
    if isinstance(conditions, list) and conditions:
        indication = "; ".join(conditions)
    else:
        indication = None

    # --- Sponsor ---
    sponsor_mod = protocol.get("sponsorCollaboratorsModule", {}) or {}
    lead = sponsor_mod.get("leadSponsor", {}) or {}
    sponsor = lead.get("name")

    # --- Results flag ---
    results_available = bool(study.get("hasResults"))

    return {
        "nct_id": nct_id,
        "phase": phase,
        "status": status,
        "indication": indication,
        "sponsor": sponsor,
        "results_available": results_available,
    }


# ------------------------------------------------------------
# Querying per ligand (by names + synonyms)
# ------------------------------------------------------------
def clinical_trials_for_names(names: list[str], max_studies_per_term: int = 200):
    """
    For all synonyms of one ligand, query ClinicalTrials.gov using query.term
    and return unique parsed study rows.
    """
    results = []
    seen_nct = set()
    searched = set()
    skipped = 0
    unrelated = 0

    for raw_name in names:
        if not isinstance(raw_name, str):
            continue
        term = raw_name.strip()
        if not term:
            continue

        # Database accessions, InChIKeys, CAS numbers and IUPAC strings are never
        # what a trial registers an intervention under, and brackets make CT.gov's
        # query parser return a 400 outright. Development codes (CD5789,
        # BMS-189961) are kept -- trials really do use those.
        term = strip_source_tag(term)
        if not is_searchable_term(term):
            skipped += 1
            continue

        key = term.lower()
        if key in searched:
            continue
        searched.add(key)

        # Pass the raw term: requests percent-encodes params itself. Encoding it
        # here first produced double-encoded queries (%27 -> %2527) that CT.gov
        # rejects with a 400, silently dropping every synonym with punctuation.
        print(f"Searching ClinicalTrials.gov for: '{term}'")
        studies = query_clinical_trials(term, max_studies=max_studies_per_term)
        time.sleep(0.25)

        for st in studies:
            parsed = parse_clinical_study(st)
            nct = parsed.get("nct_id")
            if not nct or nct in seen_nct:
                continue

            matched = study_matches_term(st, term)
            if not matched:
                unrelated += 1
                continue

            seen_nct.add(nct)
            parsed["search_term"] = term
            parsed["matched_intervention"] = matched
            results.append(parsed)

    if skipped:
        print(f"  (skipped {skipped} non-name synonym(s))")
    if unrelated:
        print(f"  (dropped {unrelated} study hit(s) with no matching intervention)")

    return results


# ------------------------------------------------------------
# Main function for your pipeline (compound_info.tsv → clinical_trials.tsv)
# ------------------------------------------------------------
def build_clinical_trials(
    compound_info_tsv: str,
    output_tsv: str = "data/processed/clinical_trials.tsv"
):
    """
    Reads data/processed/compound_info.tsv and queries CT.gov
    for name/synonym matches for each ligand.
    """
    Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(compound_info_tsv, sep="\t", dtype=str)
    rows = []

    for _, row in df.iterrows():
        ligand_id = row.get("molecule_chembl_id")
        name_field = row.get("all_pubchem_synonyms") or ""

        # Split synonyms: "A|B|C"
        synonyms = [n.strip() for n in name_field.split("|") if n.strip()]

        # Add best name first if available
        # best_name is legitimately blank for compounds with no common name, and
        # pandas gives that back as NaN, which is truthy.
        best_name = row.get("best_name") or row.get("pref_name")
        if isinstance(best_name, str) and best_name.strip():
            synonyms.insert(0, best_name.strip())

        # Deduplicate
        synonyms = list(dict.fromkeys(synonyms))

        ct_hits = clinical_trials_for_names(synonyms)

        for hit in ct_hits:
            rows.append({
                "ligand_id": ligand_id,
                "search_term": hit["search_term"],
                "nct_id": hit["nct_id"],
                "phase": hit["phase"],
                "status": hit["status"],
                "indication": hit["indication"],
                "sponsor": hit["sponsor"],
                "results_available": hit["results_available"],
                "matched_intervention": hit["matched_intervention"],
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(output_tsv, sep="\t", index=False)
    print(f"Wrote clinical trials: {df_out.shape} → {output_tsv}")
    return df_out


# ------------------------------------------------------------
# Script entry point
# ------------------------------------------------------------
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python fetch_clinical_trials.py data/processed/compound_info.tsv")
        raise SystemExit(1)

    input_path = sys.argv[1]
    build_clinical_trials(input_path)
