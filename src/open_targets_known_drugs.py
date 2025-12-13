import requests
import pandas as pd
from pathlib import Path

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

KNOWN_DRUGS_QUERY = """
query KnownDrugsForTarget($ensemblId: String!, $size: Int) {
  target(ensemblId: $ensemblId) {
    knownDrugs(size: $size) {
      count
      rows {
        phase
        status
        urls {
          name
          url
        }
        disease {
          id
          name
        }
        drug {
          id
          name
          isApproved
        }
      }
    }
  }
}
"""

def fetch_known_drugs_for_target(
    ensembl_id: str,
    size: int = 100,
    timeout: float = 30.0
) -> list[dict]:
    """
    Query Open Targets for known drugs of a target.
    Returns a list of Drug rows.
    """
    payload = {
        "query": KNOWN_DRUGS_QUERY,
        "variables": {"ensemblId": ensembl_id, "size": size}
    }
    headers = {"Content-Type": "application/json"}

    r = requests.post(API_URL, json=payload, headers=headers, timeout=timeout)
    r.raise_for_status()
    data = r.json()

    if "errors" in data:
        raise RuntimeError(f"GraphQL errors: {data['errors']}")

    tg = data.get("data", {}).get("target")
    if not tg or "knownDrugs" not in tg:
        return []
    return tg["knownDrugs"].get("rows", [])

def parse_known_drugs(raw_rows: list[dict]) -> pd.DataFrame:
    """
    Normalize raw known drugs rows into a DataFrame.
    """
    records = []
    for row in raw_rows:
        drug = row.get("drug") or {}
        disease = row.get("disease") or {}
        urls = row.get("urls") or []

        url_strs = [f"{u.get('name')}: {u.get('url')}" for u in urls if u.get("url")]

        records.append({
            "drug_id": drug.get("id"),
            "drug_name": drug.get("name"),
            "is_approved": drug.get("isApproved"),
            "phase": row.get("phase"),
            "status": row.get("status"),
            "indication_id": disease.get("id"),
            "indication_name": disease.get("name"),
            "urls": "; ".join(url_strs) if url_strs else None,
        })
    return pd.DataFrame(records)

def get_known_drugs_for_target(
    ensembl_id: str,
    size: int = 100
) -> pd.DataFrame:
    """
    Fetch and parse known drugs for a target into a pandas DataFrame.
    """
    raw = fetch_known_drugs_for_target(ensembl_id, size=size)
    return parse_known_drugs(raw)

def save_known_drugs(
    df: pd.DataFrame,
    path: str = "data/processed/open_targets_known_drugs.tsv"
) -> pd.DataFrame:
    """
    Save the known drugs DataFrame to a TSV.
    """
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    print(f"Wrote {df.shape[0]} known drugs → {path}")
    return df
