#!/usr/bin/env python3
import requests
import pandas as pd
from pathlib import Path
import time

CT_BASE = "https://clinicaltrials.gov/api/v2/studies"


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

    params = {
        "query.term": term,
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
    skipped = 0

    for raw_name in names:
        term = raw_name.strip()
        if not term:
            continue

        # CT.gov's query parser rejects brackets outright, and no trial is
        # registered under a full IUPAC name anyway. Skip instead of spending a
        # request on a guaranteed 400.
        if any(ch in term for ch in "[]"):
            skipped += 1
            continue

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

            seen_nct.add(nct)
            parsed["search_term"] = term
            results.append(parsed)

    if skipped:
        print(f"  (skipped {skipped} structural/IUPAC synonym(s) CT.gov cannot parse)")

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
        best_name = row.get("best_pubchem_name") or row.get("pref_name")
        if best_name:
            synonyms.insert(0, best_name)

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
