#!/usr/bin/env python3
import requests
import pandas as pd
from pathlib import Path
from urllib.parse import urlparse
import time
from requests.exceptions import ReadTimeout, ConnectionError, HTTPError

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PUBCHEM_PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

def parse_chembl_target_id(x: str) -> str:
    s = x.strip()
    if s.upper().startswith("CHEMBL"):
        return s.upper()
    p = urlparse(s)
    parts = p.path.rstrip("/").split("/")
    if parts:
        cand = parts[-1]
        if cand.upper().startswith("CHEMBL"):
            return cand.upper()
    raise ValueError(f"Cannot parse ChEMBL target ID from '{x}'")

def robust_get(url, params=None, headers=None,
               max_retries=5,
               timeout=(10, 120),
               allow_500_retry=True):
    """
    GET with retries on:
      - ReadTimeout
      - ConnectionError
      - HTTP 5xx (if allow_500_retry=True)

    Returns:
        requests.Response

    Raises:
        Exception after max retries
    """
    headers = headers or {}

    for attempt in range(max_retries):
        try:
            response = requests.get(
                url,
                params=params,
                headers=headers,
                timeout=timeout,
            )

            # Retry on 5xx errors
            if allow_500_retry and 500 <= response.status_code < 600:
                print(f"[WARN] Server error {response.status_code} on {url}")
                raise HTTPError(f"Server error {response.status_code}", response=response)

            # Successful response
            response.raise_for_status()
            return response

        except (ReadTimeout, ConnectionError, HTTPError) as e:
            wait = 2 ** attempt
            print(f"[Retry] Error on attempt {attempt+1}/{max_retries}: {e}. "
                  f"Retrying in {wait}s...")
            time.sleep(wait)

    # Out of retries → raise real error
    raise RuntimeError(f"Failed after {max_retries} retries → {url}")


def fetch_activities(target_chembl_id: str):
    url = f"{CHEMBL_BASE}/activity.json"
    params = {"target_chembl_id": target_chembl_id, "limit": 1000}
    activities = []
    while True:
        r = robust_get(url, params=params, timeout=(10, 120))
        j = r.json()

        batch = j.get("activities", [])
        activities.extend(batch)

        next_p = j.get("page_meta", {}).get("next")
        if not next_p:
            break

        # ChEMBL gives only the path, so prepend domain
        url = "https://www.ebi.ac.uk" + next_p
        params = None  # pagination URL already includes params
        time.sleep(0.1)
    return activities

def fetch_molecule(mol_chembl_id: str):
    url = f"{CHEMBL_BASE}/molecule/{mol_chembl_id}.json"
    r = robust_get(url, timeout=(10, 120))
    return r.json()

def pubchem_cid_from_inchikey(inchi_key: str):
    """
    Use PubChem PUG-REST to get CID for given InChIKey.
    Returns CID (int) or None if not found / error.
    """
    url = f"{PUBCHEM_PUG}/inchikey/{inchi_key}/cids/JSON"
    r = requests.get(url, timeout=120)
    if r.status_code != 200:
        return None
    j = r.json()
    cids = j.get("IdentifierList", {}).get("CID", [])
    return cids[0] if cids else None

def fetch_pubchem_metadata(cid: int):
    """
    Given a PubChem CID, fetch compound summary via PUG-REST.
    Returns dict with e.g. synonyms and other identifiers (if available).
    """
    url = f"{PUBCHEM_PUG}/cid/{cid}/record/JSON"
    r = requests.get(url, timeout=120)
    r.raise_for_status()
    return r.json()

def build_datasets(target_inputs, out_dir="data/raw"):
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    act_rows = []
    lig_rows = {}
    for t in target_inputs:
        tid = parse_chembl_target_id(t)
        print("Fetching activities for target:", tid)
        acts = fetch_activities(tid)
        print(" → got", len(acts), "activity records")

        for a in acts:
            act_rows.append({
                "target_chembl_id": tid,
                "molecule_chembl_id": a.get("molecule_chembl_id"),
                "assay_chembl_id": a.get("assay_chembl_id"),
                "standard_type": a.get("standard_type"),
                "standard_relation": a.get("standard_relation"),
                "standard_value": a.get("standard_value"),
                "standard_units": a.get("standard_units"),
                "pchembl_value": a.get("pchembl_value"),
                "activity_comment": a.get("activity_comment"),
            })

        mol_ids = {a.get("molecule_chembl_id") for a in acts if a.get("molecule_chembl_id")}
        for mol in mol_ids:
            if mol in lig_rows:
                continue
            try:
                m = fetch_molecule(mol)
            except Exception as e:
                print("Failed fetch molecule", mol, e)
                continue

            mol_struct = m.get("molecule_structures") or {}
            smiles = mol_struct.get("canonical_smiles")
            inchi_key = mol_struct.get("standard_inchi_key")

            names = set()
            pref = m.get("pref_name")
            if pref:
                names.add(pref)
            for syn in m.get("molecule_synonyms", []):
                nm = syn.get("synonym")
                if nm:
                    names.add(nm)

            pubchem_cid = None
            cas_number = None
            if inchi_key:
                cid = pubchem_cid_from_inchikey(inchi_key)
                if cid:
                    pubchem_cid = cid
                    try:
                        pc = fetch_pubchem_metadata(cid)
                        # attempt to extract CAS if present in PC data
                        # Note: PubChem does not always provide CAS RN; you may need
                        # to search synonyms or external registry IDs.
                        # Example retrieval logic (may need adaptation):
                        # pc["PC_Compounds"] → list → check each record for "RN" or "RegistryNumber"
                        # For now we skip or store None.
                    except Exception as e:
                        print("Failed fetch PubChem metadata for CID", cid, e)
                time.sleep(0.35)

            lig_rows[mol] = {
                "molecule_chembl_id": mol,
                "canonical_smiles": smiles,
                "standard_inchi_key": inchi_key,
                "names": "|".join(sorted(names)),
                "pubchem_cid": pubchem_cid,
                "cas_number": cas_number,
            }

        time.sleep(0.1)

    df_act = pd.DataFrame(act_rows)
    df_lig = pd.DataFrame(list(lig_rows.values()))

    df_act.to_csv(Path(out_dir) / "activities.tsv", sep="\t", index=False)
    df_lig.to_csv(Path(out_dir) / "ligands.tsv", sep="\t", index=False)

    print("Wrote activities:", df_act.shape, "ligands:", df_lig.shape)
    return df_act, df_lig
