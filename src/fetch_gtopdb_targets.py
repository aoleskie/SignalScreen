#!/usr/bin/env python3
# fetch_gtopdb_target.py

import re
import requests
from pathlib import Path
import pandas as pd
from bs4 import BeautifulSoup
import time

API_BASE = "https://www.guidetopharmacology.org/services/targets"
WEB_BASE = "https://www.guidetopharmacology.org/GRAC/ObjectDisplayForward"

def parse_gtop_id(x: str) -> int:
    x = x.strip()
    # bare numeric
    m = re.fullmatch(r"\d+", x)
    if m:
        return int(x)
    # objectId= in URL
    m2 = re.search(r"[?&]objectId=(\d+)", x)
    if m2:
        return int(m2.group(1))
    raise ValueError(f"Cannot parse GtoPdb target ID from '{x}'")

def single_target_api(tid: int) -> dict:
    url = f"{API_BASE}/{tid}"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    return resp.json()

def synonym_target_api(tid: int) -> dict:
    url = f"{API_BASE}/{tid}/synonyms"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    return resp.json()

def genename_target_api(tid: int) -> dict:
    url = f"{API_BASE}/{tid}/geneProteinInformation?species=Human"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    return resp.json()

def database_target_api(tid: int) -> dict:
    url = f"{API_BASE}/{tid}/databaseLinks?species=Human"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    return resp.json()

def extract_db_links(db_links_json, want_databases=None):
    """
    From a list of database-link dicts (as from GtoPdb /databaseLinks),
    return a dict mapping each database name in want_databases to its URL (or accession if URL missing).
    If want_databases is None, return a dict of all database → url/accession.
    If no links or no matching entries, return {}.
    """
    if not db_links_json:
        return {}

    result = {}
    for entry in db_links_json:
        db = entry.get("database")
        # normalize database name if needed
        url = entry.get("url") or entry.get("accession")
        if not db:
            continue
        if want_databases:
            if db in want_databases:
                result[db] = url
        else:
            result[db] = url
    return result

def scrape_target_page(tid: int) -> dict:
    url = f"{WEB_BASE}?objectId={tid}"
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    soup = BeautifulSoup(resp.text, "html.parser")
    out = {}
    # Systematic Nomenclature
    b_tag = soup.find("b", string=re.compile(r"\bSystematic Nomenclature\b"))
    if b_tag and b_tag.next_sibling:
        out["systematic_nomenclature"] = b_tag.next_sibling.strip() or None
    else:
        out["systematic_nomenclature"] = None

    return out

def fetch_target_full(tid: int) -> dict:
    #get single target api, name, type, familyIDs
    rec = {"gtop_target_id": tid}
    try:
        api = single_target_api(tid)
    except Exception as e:
        print(f"API fetch failed for {tid}: {e}")
        api = {}
    # copy basic fields
    rec["name"] = api.get("name")
    rec["type"] = api.get("type")
    rec["familyIds"] = "|".join(str(i) for i in api.get("familyIds", [])) if api.get("familyIds") else None

    #get synonyms
    try:
        api = synonym_target_api(tid)
    except Exception as e:
        print(f"Synonym fetch failed for {tid}: {e}")
        api = {}
    if not api:
        rec["synonyms"] = None
    else:
        rec["synonyms"] = [d.get("name") for d in api if d.get("name") is not None]
        # print(rec["synonyms"])

    #get gene name
    try:
        api = genename_target_api(tid)[0]
    except Exception as e:
        print(f"Gene symbol fetch failed for {tid}: {e}")
        api = {}
    if not api:
        rec["gene_symbol"] = None
    else:
        rec["gene_symbol"] = api.get("geneSymbol")

    #get DataBase links
    try:
        api = database_target_api(tid)
    except Exception as e:
        print(f"Database fetch failed for {tid}: {e}")
        api = {}
    if api:
        db_map = extract_db_links(api, want_databases=["ChEMBL Target", "DrugBank Target"])
        rec["chembl_target_url"] = db_map.get("ChEMBL Target")
        rec["drugbank_target_url"] = db_map.get("DrugBank Target")



    # fallback data from HTML
    try:
        html = scrape_target_page(tid)
    except Exception as e:
        print(f"HTML scrape failed for {tid}: {e}")
        html = {}

    # Merge: prefer API-derived or HTML fallback if missing
    rec["systematic_nomenclature"] = html.get("systematic_nomenclature")
    return rec

def fetch_targets(ids_or_urls, output_basename="gtop_targets_full.tsv"):
    rows = []
    for x in ids_or_urls:
        try:
            tid = parse_gtop_id(x)
        except ValueError as e:
            print("Warning:", e)
            continue
        print("Processing target:", tid)
        rec = fetch_target_full(tid)
        rows.append(rec)
        # polite delay to avoid hammering server
        time.sleep(0.5)

    if not rows:
        print("No targets processed.")
        return None

    df = pd.DataFrame(rows)
    out_dir = Path.cwd() / "data" / "targets"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / output_basename
    df.to_csv(out_path, sep="\t", index=False)
    print("Wrote", len(df), "entries to", out_path)
    return df

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fetch GtoPdb target metadata (API + HTML fallback)")
    parser.add_argument("inputs", nargs="+", help="Target ID(s) or URL(s)")
    parser.add_argument("--out", "-o", default="gtop_targets_full.tsv", help="Output TSV basename")
    args = parser.parse_args()
    fetch_targets(args.inputs, output_basename=args.out)
