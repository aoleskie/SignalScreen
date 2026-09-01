#!/usr/bin/env python3

import requests
import pandas as pd
from pathlib import Path
import re
import time

from src.compound_names import choose_display_name, official_name, is_bare_code

PUGREST_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_BATCH_SIZE = 25


def fetch_chembl_pref_names(chembl_ids):
    """
    ChEMBL's curated pref_name for each molecule, batched.

    This is the best available name when it exists: 'TRETINOIN' rather than
    whichever synonym PubChem happens to list first.
    """
    out = {}
    ids = [i for i in chembl_ids if isinstance(i, str) and i.strip()]
    for i in range(0, len(ids), CHEMBL_BATCH_SIZE):
        batch = ids[i:i + CHEMBL_BATCH_SIZE]
        try:
            r = requests.get(
                f"{CHEMBL_BASE}/molecule.json",
                params={
                    "molecule_chembl_id__in": ",".join(batch),
                    "limit": CHEMBL_BATCH_SIZE,
                    "only": "molecule_chembl_id,pref_name",
                },
                timeout=60,
            )
            r.raise_for_status()
        except Exception as e:
            print(f"[warn] ChEMBL pref_name batch failed: {e}")
            continue
        for m in r.json().get("molecules", []):
            if m.get("pref_name"):
                out[m["molecule_chembl_id"]] = m["pref_name"]
        time.sleep(0.2)
    return out

def get_pubchem_cid_from_inchikey(inchi_key: str):
    """Map an InChIKey to PubChem CID via PUG-REST."""
    if not inchi_key or pd.isna(inchi_key):
        return None
    url = f"{PUGREST_BASE}/inchikey/{inchi_key}/cids/JSON"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
    except Exception:
        return None
    j = r.json()
    return j.get("IdentifierList", {}).get("CID", [None])[0]

def fetch_pubchem_synonyms(cid: int):
    """Retrieve synonym list for given PubChem CID."""
    if not cid:
        return []
    url = f"{PUGREST_BASE}/cid/{cid}/synonyms/JSON"
    try:
        r = requests.get(url, timeout=15)
        r.raise_for_status()
    except Exception:
        return []
    j = r.json()
    infos = j.get("InformationList", {}).get("Information", [])
    if not infos:
        return []
    return infos[0].get("Synonym", [])

def extract_cas_from_synonyms(synonyms):
    """Extract a CAS-like string from synonyms, if present."""
    cas_pattern = re.compile(r"\b\d{2,7}-\d{2}-\d\b")
    for s in synonyms:
        if cas_pattern.search(s):
            return s
    return None

def build_compound_info(
    sorted_tsv: str = "data/processed/sorted_ligands_RARG.tsv",
    raw_ligands_tsv: str = "data/raw/ligands.tsv",
    n: int = 10,
    out_dir: str = "data/processed"
):
    """
    Enrich the top N sorted ligands with PubChem CID and synonyms,
    outputting compound_info.tsv in data/processed.
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    # Load sorted (top candidates)
    df_sorted = pd.read_csv(sorted_tsv, sep="\t")
    top_ids = df_sorted["molecule_chembl_id"].tolist()[:n]

    # Load full raw ligand info
    df_raw = pd.read_csv(raw_ligands_tsv, sep="\t")

    # Join sorted -> raw on ChEMBL ID to get InChIKey and SMILES
    df_top = pd.merge(
        pd.DataFrame({"molecule_chembl_id": top_ids}),
        df_raw,
        how="left",
        on="molecule_chembl_id"
    )

    pref_names = fetch_chembl_pref_names(top_ids)

    records = []
    for _, row in df_top.iterrows():
        chembl_id = row["molecule_chembl_id"]
        smiles = row.get("canonical_smiles")
        inchi_key = row.get("standard_inchi_key")

        pubchem_cid = None
        all_synonyms = None
        cas_number = None
        syns = []

        if pd.notna(inchi_key):
            pubchem_cid = get_pubchem_cid_from_inchikey(inchi_key)
            time.sleep(0.1)

            if pubchem_cid:
                syns = fetch_pubchem_synonyms(pubchem_cid)
                time.sleep(0.1)
                if syns:
                    all_synonyms = "|".join(syns)
                    cas_number = extract_cas_from_synonyms(syns)

        # Name preference: ChEMBL's curated pref_name, then an INN/USAN-tagged
        # synonym, then the most readable remaining synonym. None of those
        # existing means the compound genuinely has no common name -- leave it
        # blank rather than presenting an accession as a name.
        best_name = pref_names.get(chembl_id)
        name_source = "chembl_pref_name" if best_name else None

        if not best_name:
            best_name = official_name(syns)
            name_source = "inn_usan" if best_name else None

        if not best_name:
            best_name = choose_display_name(syns)
            if not best_name:
                name_source = "unnamed"
            elif is_bare_code(best_name):
                # A catalogue/development code is the best available label, but
                # it is not a name; say so rather than implying otherwise.
                name_source = "code"
            else:
                name_source = "synonym"

        records.append({
            "molecule_chembl_id": chembl_id,
            "canonical_smiles": smiles,
            "standard_inchi_key": inchi_key,
            "pubchem_cid": pubchem_cid,
            "best_name": best_name,
            "name_source": name_source,
            "all_pubchem_synonyms": all_synonyms,
            "cas_number": cas_number,
        })

    out_df = pd.DataFrame(records)
    out_path = Path(out_dir) / "compound_info.tsv"
    out_df.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote {len(out_df)} compounds to {out_path}")

    return out_df
