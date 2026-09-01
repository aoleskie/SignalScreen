#!/usr/bin/env python3
import requests
import pandas as pd
import re
from pathlib import Path
import time
import urllib.parse

NCBI_EUTIL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"


# ------------------------------------------------------------
# ROBUST GET
# ------------------------------------------------------------
def robust_get(
    url,
    params=None,
    timeout=45,
    retries=4,
    backoff=1.5,
    ok_status=(200,),
    verbose=True,
):
    last_exc = None

    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, params=params, timeout=timeout)
            if r.status_code in ok_status:
                return r

            if verbose:
                print(
                    f"[warn] GET {url} → HTTP {r.status_code} "
                    f"(attempt {attempt}/{retries})"
                )

        # Broad on purpose: a truncated body (ChunkedEncodingError /
        # ProtocolError) is just as transient as a timeout, and letting it escape
        # discards the whole run's work.
        except requests.exceptions.RequestException as e:
            last_exc = e
            if verbose:
                print(
                    f"[retry] GET {url} failed ({e.__class__.__name__}) "
                    f"(attempt {attempt}/{retries})"
                )

        time.sleep(backoff * attempt)

    raise RuntimeError(f"GET failed after {retries} attempts: {url}") from last_exc


# ------------------------------------------------------------
# ASSAY VALUE EXTRACTION
# ------------------------------------------------------------
ASSAY_PATTERN = re.compile(
    r"(?P<type>IC50|EC50|Ki|Kd)\s*[=:\s]\s*(?P<value>[0-9.]+(?:e-?\d+)?)\s*(?P<unit>nM|uM|μM|pM)?",
    re.IGNORECASE
)

def extract_assay_values(text: str):
    results = []
    for m in ASSAY_PATTERN.finditer(text or ""):
        results.append({
            "assay_type": m.group("type"),
            "assay_value": m.group("value"),
            "units": m.group("unit")
        })
    return results


# ------------------------------------------------------------
# PUBMED SEARCH (FIRST-PASS NAME SEARCH ONLY)
# ------------------------------------------------------------
def pubmed_esearch_simple(term: str, retmax=200):
    """
    Simple PubMed ESearch:
    https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=PALOVAROTENE
    """
    if not term or not term.strip():
        return []

    params = {"db": "pubmed", "term": term, "retmax": retmax, "retmode": "json"}
    r = robust_get(f"{NCBI_EUTIL}/esearch.fcgi", params=params, timeout=20)
    data = r.json()
    return data.get("esearchresult", {}).get("idlist", [])


EFETCH_BATCH_SIZE = 100


def pubmed_efetch(pmids: list[str]):
    """Fetch PubMed records, batched. One 200-ID request is more likely to come
    back truncated than two 100-ID ones."""
    if not pmids:
        return []

    documents = []
    for i in range(0, len(pmids), EFETCH_BATCH_SIZE):
        batch = pmids[i:i + EFETCH_BATCH_SIZE]
        params = {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"}
        try:
            r = robust_get(f"{NCBI_EUTIL}/efetch.fcgi", params=params, timeout=60)
        except RuntimeError as e:
            print(f"[warn] efetch batch {i // EFETCH_BATCH_SIZE} failed, skipping: {e}")
            continue
        documents.append(r.text)
        time.sleep(0.34)  # NCBI allows 3 requests/sec without an API key
    return documents


def parse_pubmed_xml(xml_docs):
    import xml.etree.ElementTree as ET
    if not xml_docs:
        return []
    if isinstance(xml_docs, str):
        xml_docs = [xml_docs]

    articles = []
    for xml_text in xml_docs:
        if not xml_text:
            continue
        try:
            articles.extend(ET.fromstring(xml_text).findall(".//PubmedArticle"))
        except ET.ParseError as e:
            print(f"[warn] unparseable PubMed XML batch, skipping: {e}")

    papers = []

    for article in articles:
        pmid = article.findtext(".//PMID")
        title = article.findtext(".//ArticleTitle")

        # extract year
        year = article.findtext(".//PubDate/Year")
        if not year:
            md = article.findtext(".//PubDate/MedlineDate")
            year = md[:4] if md else None

        abstract = " ".join([elem.text or "" for elem in article.findall(".//AbstractText")])

        authors_list = []
        for a in article.findall(".//Author"):
            last = a.findtext("LastName") or ""
            fore = a.findtext("ForeName") or ""
            if last:
                authors_list.append(f"{last} {fore}".strip())
        authors = "; ".join(authors_list)

        assays = extract_assay_values(abstract)

        base = {
            "pmid": pmid,
            "title": title,
            "year": year,
            "authors": authors,
            "abstract": abstract,
            "link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
        }

        if assays:
            for a in assays:
                papers.append({
                    **base,
                    "assay_type": a["assay_type"],
                    "assay_value": a["assay_value"],
                    "units": a["units"],
                })
        else:
            papers.append({
                **base,
                "assay_type": None,
                "assay_value": None,
                "units": None,
            })

    return papers


# ------------------------------------------------------------
# BINDINGDB
# ------------------------------------------------------------
def query_bindingdb_inchikey(inchikey: str):
    if not inchikey:
        return []

    url = (
        "https://www.bindingdb.org/axis2/services/BDBRest?"
        f"operation=getLigandByInchiKey&InchiKey={inchikey}"
    )
    try:
        r = robust_get(url, timeout=30, retries=1)
        text = r.text

        rows = []
        for line in text.splitlines():
            if line.startswith("Ligand"):
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            rows.append({
                "bindingdb_id": parts[0],
                "target_protein": parts[1],
                "affinity": parts[2],
                "units": parts[3],
                "source_pmid": parts[4] if parts[4] != "NA" else None,
            })
        return rows

    except Exception as e:
        print(f"[BindingDB] error for InChIKey={inchikey}: {e}")
        return []


# ------------------------------------------------------------
# ChEMBL DOCUMENT LINKAGE
# ------------------------------------------------------------
DOCUMENT_BATCH_SIZE = 20


def fetch_activity_documents_for_molecule(mol_chembl_id: str):
    """
    Collect the ChEMBL document IDs behind a molecule's activities.

    The activity resource already carries document_chembl_id, so there is no need
    to resolve each activity's assay first. Retinoic acid has 2332 activity rows
    over 1843 assays; the per-assay walk cost ~4000 sequential requests for that
    one molecule, which this reduces to one paginated pass.
    """
    docs = set()
    url = f"{CHEMBL_API}/activity.json"
    params = {
        "molecule_chembl_id": mol_chembl_id,
        "limit": 1000,
        "only": "document_chembl_id",
    }

    while True:
        r = robust_get(url, params=params, timeout=120)
        data = r.json()

        for act in data.get("activities", []):
            doc_id = act.get("document_chembl_id")
            if doc_id:
                docs.add(doc_id)

        next_page = data.get("page_meta", {}).get("next")
        if not next_page:
            break

        url = "https://www.ebi.ac.uk" + next_page
        params = None

    return sorted(docs)


def fetch_chembl_document_records(document_ids: list[str]):
    """Fetch document metadata in batches via document_chembl_id__in."""
    out = []
    for i in range(0, len(document_ids), DOCUMENT_BATCH_SIZE):
        batch = document_ids[i:i + DOCUMENT_BATCH_SIZE]
        r = robust_get(
            f"{CHEMBL_API}/document.json",
            params={
                "document_chembl_id__in": ",".join(batch),
                "limit": DOCUMENT_BATCH_SIZE,
            },
            timeout=120,
        )
        for doc in r.json().get("documents", []):
            doc_id = doc.get("document_chembl_id")
            pmid = doc.get("pubmed_id")
            out.append({
                "chembl_document_id": doc_id,
                "pmid": str(pmid) if pmid else None,
                "title": doc.get("title"),
                "year": doc.get("year"),
                "doi": doc.get("doi"),
                "chembl_doc_link": f"https://www.ebi.ac.uk/chembl/document_report_card/{doc_id}/",
            })
    return out


# ------------------------------------------------------------
# MAIN: unified literature builder
# ------------------------------------------------------------
def build_literature(
    compound_info_tsv="data/processed/compound_info.tsv",
    ligand_targets_tsv="data/processed/ligand_targets.tsv",
    out_merged="data/processed/literature_merged.tsv",
):
    lig = pd.read_csv(compound_info_tsv, sep="\t", dtype=str)

    rows = []

    for _, row in lig.iterrows():
        ligand_id = row["molecule_chembl_id"]
        print(f"Literature lookup for {ligand_id} ...")

        # Collect synonyms
        raw_names = row.get("all_pubchem_synonyms") or ""
        name_list = [n.strip() for n in raw_names.split("|") if n.strip()]

        for extra in (row.get("best_pubchem_name"), row.get("pref_name")):
            if extra and isinstance(extra, str) and extra.strip():
                name_list.append(extra.strip())

        name_list = list(dict.fromkeys(name_list))  # dedupe

        # ------------------------------------------------------------
        # PUBMED FIRST-PASS SEARCH (term = NAME)
        # ------------------------------------------------------------
        seen_pmids = set()

        for term in name_list:
            print(term)
            try:
                pmids = pubmed_esearch_simple(term)
                if not pmids:
                    continue
                papers = parse_pubmed_xml(pubmed_efetch(pmids))
            except Exception as e:
                print(f"[warn] PubMed lookup failed for {ligand_id} / '{term}': {e}")
                continue

            for p in papers:
                pmid = p["pmid"]
                if pmid in seen_pmids:
                    continue
                seen_pmids.add(pmid)

                rows.append({
                    "ligand_id": ligand_id,
                    "pmid": pmid,
                    "title": p["title"],
                    "year": p["year"],
                    "assay_type": p.get("assay_type"),
                    "assay_value": p.get("assay_value"),
                    "units": p.get("units"),
                    "abstract": p.get("abstract"),
                    "target": None,
                    "affinity": None,
                    "doi": None,
                    "source": "PubMed",
                    "link": p.get("link"),
                })

        # ------------------------------------------------------------
        # BINDINGDB
        # ------------------------------------------------------------
        inchikey = row.get("standard_inchi_key")
        if inchikey:
            bdb_hits = query_bindingdb_inchikey(inchikey)
            for b in bdb_hits:
                rows.append({
                    "ligand_id": ligand_id,
                    "pmid": b.get("source_pmid"),
                    "title": None,
                    "year": None,
                    "assay_type": None,
                    "assay_value": b.get("affinity"),
                    "units": b.get("units"),
                    "abstract": None,
                    "target": b.get("target_protein"),
                    "affinity": b.get("affinity"),
                    "doi": None,
                    "source": "BindingDB",
                    "link": (
                        f"https://pubmed.ncbi.nlm.nih.gov/{b['source_pmid']}/"
                        if b.get("source_pmid") else None
                    ),
                })

        # ------------------------------------------------------------
        # ChEMBL LITERATURE
        # ------------------------------------------------------------
        try:
            doc_ids = fetch_activity_documents_for_molecule(ligand_id)
            docs = fetch_chembl_document_records(doc_ids) if doc_ids else []
        except Exception as e:
            print(f"[warn] ChEMBL literature failed for {ligand_id}: {e}")
            docs = []

        if docs:
            for c in docs:
                rows.append({
                    "ligand_id": ligand_id,
                    "pmid": c.get("pmid"),
                    "title": c.get("title"),
                    "year": c.get("year"),
                    "assay_type": None,
                    "assay_value": None,
                    "units": None,
                    "abstract": None,
                    "target": None,
                    "affinity": None,
                    "doi": c.get("doi"),
                    "source": "ChEMBL",
                    "link": (
                        f"https://pubmed.ncbi.nlm.nih.gov/{c['pmid']}/"
                        if c.get("pmid") else c.get("chembl_doc_link")
                    ),
                })

        time.sleep(0.25)

    df_out = pd.DataFrame(rows)
    Path(out_merged).parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(out_merged, sep="\t", index=False)

    print(f"Wrote unified literature table: {df_out.shape} → {out_merged}")
    return df_out
