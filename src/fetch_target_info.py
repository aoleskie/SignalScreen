import requests
import pandas as pd
from pathlib import Path
import time

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
REACTOME_BASE = "https://reactome.org/AnalysisService/identifiers"

# Only these assay endpoints are comparable as "how tightly does this bind".
# Everything else ChEMBL returns for a molecule (ADMET readouts such as Fu, Ratio
# EC50, Efficacy, activity percentages) lives on an unrelated scale and must never
# compete for "best affinity" -- a fraction-unbound of 0.0001 is not a 0.1 pM binder.
AFFINITY_TYPES = {"IC50", "EC50", "Ki", "Kd", "AC50", "XC50"}

# Multipliers onto nM, so values from different assays are actually comparable.
UNIT_TO_NM = {
    "M": 1e9,
    "mM": 1e6,
    "uM": 1e3,
    "µM": 1e3,
    "μM": 1e3,
    "nM": 1.0,
    "pM": 1e-3,
    "fM": 1e-6,
}

# A ">" relation means the assay never reached an endpoint, so the true value is
# unbounded above; treating it as an exact measurement invents potency.
EXACT_RELATIONS = {"=", None, ""}


def normalize_affinity_nm(act):
    """
    Return the activity's value in nM, or None if it is not a comparable
    affinity measurement.
    """
    if act.get("standard_type") not in AFFINITY_TYPES:
        return None
    if act.get("standard_relation") not in EXACT_RELATIONS:
        return None

    factor = UNIT_TO_NM.get(act.get("standard_units"))
    if factor is None:
        return None

    try:
        value = float(act.get("standard_value"))
    except (TypeError, ValueError):
        return None
    if value <= 0:
        return None

    return value * factor

def robust_get(
    url,
    params=None,
    timeout=40,
    retries=6,
    backoff=1.5,
    allowed_status=(200,)
):
    """
    Robust GET with retries and exponential backoff.

    - timeout: per-request timeout (seconds)
    - retries: total attempts
    - backoff: multiplier for sleep between retries
    """
    last_exc = None

    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code in allowed_status:
                return r
            else:
                print(
                    f"[warn] GET {url} returned {r.status_code} "
                    f"(attempt {attempt}/{retries})"
                )
        except (requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectTimeout,
                requests.exceptions.ConnectionError) as e:
            last_exc = e
            print(
                f"[retry] GET {url} failed ({e.__class__.__name__}) "
                f"(attempt {attempt}/{retries})"
            )

        time.sleep(backoff * attempt)

    raise RuntimeError(f"GET failed after {retries} attempts: {url}") from last_exc

def fetch_activity_for_ligand(mol_id):
    url = f"{CHEMBL_BASE}/activity.json"
    params = {"molecule_chembl_id": mol_id, "limit": 1000}
    activities = []
    while True:
        r = robust_get(url, params=params, timeout=20)
        data = r.json()
        activities.extend(data.get("activities", []))
        next_page = data.get("page_meta", {}).get("next")
        if not next_page:
            break
        url = "https://www.ebi.ac.uk" + next_page
        params = None
        time.sleep(0.1)
    return activities

def fetch_mechanisms_for_ligand(mol_id):
    url = f"{CHEMBL_BASE}/mechanism.json"
    params = {"molecule_chembl_id": mol_id}
    r = robust_get(url, params=params, timeout=20)
    return r.json().get("mechanisms", [])

def fetch_target_metadata(target_chembl_id):
    url = f"{CHEMBL_BASE}/target/{target_chembl_id}.json"
    r = robust_get(url, timeout=20)
    if r.status_code != 200:
        return {}
    data = r.json()
    return data

def map_uniprot_to_hgnc(uniprot_id):
    """Map UniProt accession to HGNC symbol via UniProt mapping (or HGNC API)."""
    if not uniprot_id:
        return None
    # UniProt ID mapping service
    # Example: https://rest.uniprot.org/idmapping/run
    # For now, return the UniProt ID as placeholder if mapping fails.
    return uniprot_id

def get_ligand_targets(mol_id):
    activities = fetch_activity_for_ligand(mol_id)
    mechanisms = fetch_mechanisms_for_ligand(mol_id)

    mech_map = {m.get("target_chembl_id"): m.get("action_type") for m in mechanisms}

    # Collect one row per target
    target_map = {}
    for act in activities:
        tid = act.get("target_chembl_id")
        if not tid:
            continue
        rec = target_map.setdefault(tid, {
            "ligand_id": mol_id,
            "target_chembl_id": tid,
            "target_name": None,
            "gene_symbol": None,
            "moa": mech_map.get(tid),
            "best_assay_type": None,
            "best_affinity": None,
            "units": None,
            "n_affinity_measurements": 0,
        })

        # Choose best potency: strongest (lowest) value among comparable
        # affinity endpoints only, after converting everything to nM.
        value_nm = normalize_affinity_nm(act)
        if value_nm is None:
            continue

        rec["n_affinity_measurements"] += 1
        current = rec["best_affinity"]
        if current is None or value_nm < current:
            rec.update({
                "best_assay_type": act.get("standard_type"),
                "best_affinity": value_nm,
                "units": "nM",
            })

    # Fetch metadata for each target
    for tid, rec in target_map.items():
        meta = fetch_target_metadata(tid)
        rec["target_name"] = meta.get("pref_name")
        # Try to get UniProt accession
        comps = meta.get("target_components", [])
        gene_sym = None
        for comp in comps:
            acc = comp.get("accession")
            if acc:
                gene_sym = map_uniprot_to_hgnc(acc)
                if gene_sym:
                    break
        rec["gene_symbol"] = gene_sym

    return list(target_map.values())

def build_target_info_from_compounds(compound_info_tsv, out_dir="data/processed"):
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    lig_df = pd.read_csv(compound_info_tsv, sep="\t")
    all_targets = []
    for _, row in lig_df.iterrows():
        mol = row["molecule_chembl_id"]
        print("Fetching targets for", mol)
        tgt_rows = get_ligand_targets(mol)
        all_targets.extend(tgt_rows)
        time.sleep(0.2)

    df_targets = pd.DataFrame(all_targets)
    out_file = Path(out_dir) / "ligand_targets.tsv"
    df_targets.to_csv(out_file, sep="\t", index=False)
    return df_targets


def fetch_reactome_pathways_for_identifiers(identifiers):
    """
    Submit a list of identifiers (UniProt or gene symbols) to the Reactome
    AnalysisService identifiers endpoint and return pathway associations.
    """

    # The API expects plain text body with one identifier per line
    body = "\n".join(identifiers)
    headers = {"Content-Type": "text/plain"}

    # Submit identifiers → receive summary token and partial result
    r = requests.post(f"{REACTOME_BASE}/projection", headers=headers, data=body, timeout=120)
    r.raise_for_status()
    result = r.json()

    # "pathways" in result contains hit pathways matching submitted IDs
    pathways = []
    for p in result.get("pathways", []):
        pathways.append({
            "pathway_id": p.get("stId"),
            "pathway_name": p.get("name"),
            "entities_found": p.get("entities", {}).get("total"),
        })
    return pathways


def build_pathway_info_from_targets(targets_tsv, out_dir="data/processed"):
    """
    For each ligand target (with gene_symbol or UniProt), call Reactome's
    AnalysisService /identifiers to map proteins to pathways.
    """
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    tgt_df = pd.read_csv(targets_tsv, sep="\t")
    pathways_rows = []

    # Reactome wants a flat list of identifiers
    # But we will group them by ligand
    grouped = tgt_df.groupby("ligand_id")

    for ligand_id, subdf in grouped:
        # Collect unique identifiers (prefer UniProt, else gene_symbol)
        ids = subdf["gene_symbol"].dropna().unique().tolist()
        if not ids:
            continue

        print(f"Querying Reactome for ligand {ligand_id} identifiers: {ids}")
        try:
            pws = fetch_reactome_pathways_for_identifiers(ids)
        except Exception as e:
            print(f"Reactome API failed for {ligand_id}: {e}")
            continue

        for pw in pws:
            pathways_rows.append({
                "ligand_id": ligand_id,
                "pathway_id": pw["pathway_id"],
                "pathway_name": pw["pathway_name"],
                "entities_found": pw["entities_found"],
            })

        # Be polite
        time.sleep(0.2)

    df_pathways = pd.DataFrame(pathways_rows)
    out_file = Path(out_dir) / "ligand_pathways.tsv"
    df_pathways.to_csv(out_file, sep="\t", index=False)
    return df_pathways
