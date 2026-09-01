#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
from pandas.errors import EmptyDataError

RARG_CHEMBL_ID = "CHEMBL2003"


# --------------------------------------------------------------------
# Helper summaries
# --------------------------------------------------------------------
def summarize_publications(lit_df, max_items=15):
    """
    Return semicolon-separated literature entries across all sources.

    Format:
      Title [Source] (link_or_None)

    Behavior:
    - Includes PubMed, ChEMBL, BindingDB, etc.
    - Keeps titles even if no link exists
    - Constructs PubMed links if PMID exists
    - Deduplicates by PMID if present, else by (title, source)
    """
    if lit_df is None or lit_df.empty:
        return None

    df = lit_df.copy()

    if "title" not in df.columns:
        return None

    # Normalize text fields
    df["title"] = df["title"].astype(str).str.strip()
    df = df[df["title"].notna() & (df["title"] != "")]

    # Ensure required columns exist
    for col in ["source", "link", "pmid"]:
        if col not in df.columns:
            df[col] = None

    # Construct PubMed links where possible
    pmid_series = df["pmid"].astype(str).str.strip()
    has_pmid = pmid_series.notna() & ~pmid_series.isin(["", "None", "nan"])
    missing_link = df["link"].isna() | (df["link"].astype(str).str.strip() == "")
    df.loc[missing_link & has_pmid, "link"] = (
        "https://pubmed.ncbi.nlm.nih.gov/" + pmid_series[missing_link & has_pmid] + "/"
    )

    # Deduplication strategy
    # 1) PMID if present
    # 2) else (title, source)
    df["pmid_norm"] = pmid_series.where(has_pmid, None)

    with_pmid = df[df["pmid_norm"].notna()].drop_duplicates(subset=["pmid_norm"])
    no_pmid = df[df["pmid_norm"].isna()].drop_duplicates(subset=["title", "source"])
    df = pd.concat([with_pmid, no_pmid], ignore_index=True)

    # Prefer entries with links
    df["_has_link"] = ~(df["link"].isna() | (df["link"].astype(str).str.strip() == ""))
    df = df.sort_values("_has_link", ascending=False)

    # Format output
    items = []
    for _, row in df.head(max_items).iterrows():
        title = row["title"]
        source = row.get("source") or "Unknown"
        link = row.get("link")
        link_str = link if (link and str(link).strip() != "") else None
        items.append(f"{title} [{source}] ({link_str})")

    return "; ".join(items) if items else None




def summarize_clinical_trials(clin_df, max_items=3):
    """
    Return semicolon-separated 'NCTxxxx (clinicaltrials.gov link)' strings.
    """
    if clin_df is None or clin_df.empty:
        return None

    trials = (
        clin_df.dropna(subset=["nct_id"])
        .drop_duplicates(subset=["nct_id"])
    )

    items = [
        f"{row['nct_id']} (https://clinicaltrials.gov/study/{row['nct_id']})"
        for _, row in trials.head(max_items).iterrows()
    ]

    return "; ".join(items) if items else None


# --------------------------------------------------------------------
# Main summary builder
# --------------------------------------------------------------------
def build_ligand_summary(
    compound_info="data/processed/compound_info.tsv",
    ligand_targets="data/processed/ligand_targets.tsv",
    ligand_pathways="data/processed/ligand_pathways.tsv",
    clinical_trials="data/processed/clinical_trials.tsv",
    literature="data/processed/literature_merged.tsv",
    output="data/processed/ligand_deep_dive_summary.tsv",
):
    # ----------------------------------------------------------------
    # Load required tables
    # ----------------------------------------------------------------
    df_compound = pd.read_csv(compound_info, sep="\t", dtype=str)
    df_targets = pd.read_csv(ligand_targets, sep="\t", dtype=str)
    df_paths = pd.read_csv(ligand_pathways, sep="\t", dtype=str)
    df_lit = pd.read_csv(literature, sep="\t", dtype=str)

    # --- Clinical trials: robust to empty or missing file ---
    try:
        df_clin = pd.read_csv(clinical_trials, sep="\t", dtype=str)
        if df_clin.empty:
            raise EmptyDataError("clinical_trials.tsv is empty")
    except (FileNotFoundError, EmptyDataError) as e:
        print(f"[INFO] No clinical trials found ({e}). Using empty table.")
        df_clin = pd.DataFrame(columns=["ligand_id", "nct_id"])

    # ----------------------------------------------------------------
    # Audit
    # ----------------------------------------------------------------
    print("\n=== Column Audit ===")
    print("compound_info columns:", list(df_compound.columns))
    print("ligand_targets columns:", list(df_targets.columns))
    print("ligand_pathways columns:", list(df_paths.columns))
    print("clinical_trials columns:", list(df_clin.columns))
    print("literature_merged columns:", list(df_lit.columns))
    print("====================\n")

    # ----------------------------------------------------------------
    # Pre-index groupings
    # ----------------------------------------------------------------
    targets_by_lig = df_targets.groupby("ligand_id")
    paths_by_lig = df_paths.groupby("ligand_id")
    lit_by_lig = df_lit.groupby("ligand_id")
    clin_by_lig = df_clin.groupby("ligand_id") if not df_clin.empty else {}

    summary_rows = []

    # ----------------------------------------------------------------
    # Helper: best target per ligand (lowest affinity)
    # ----------------------------------------------------------------
    def best_target_for_lig(lig):
        """Strongest-binding target for a ligand. Affinities are nM from
        comparable endpoints only (see fetch_target_info.normalize_affinity_nm),
        so the minimum is meaningful across assays."""
        if lig not in targets_by_lig.groups:
            return None, None, None

        sub = targets_by_lig.get_group(lig).copy()
        if "best_affinity" not in sub.columns:
            return None, None, None

        sub["numeric_aff"] = pd.to_numeric(sub["best_affinity"], errors="coerce")
        sub = sub.dropna(subset=["numeric_aff"])
        if sub.empty:
            return None, None, None

        top = sub.sort_values("numeric_aff", ascending=True).iloc[0]
        gene = top.get("gene_symbol")
        name = top.get("target_name")

        return f"{gene} ({name})", top["numeric_aff"], top.get("best_assay_type")

    def rarg_affinity_for_lig(lig):
        """RARG affinity in nM, which is what the old best_rarg_affinity column
        claimed to hold but did not."""
        if lig not in targets_by_lig.groups:
            return None, None
        sub = targets_by_lig.get_group(lig)
        sub = sub[sub["target_chembl_id"] == RARG_CHEMBL_ID].copy()
        if sub.empty or "best_affinity" not in sub.columns:
            return None, None
        sub["numeric_aff"] = pd.to_numeric(sub["best_affinity"], errors="coerce")
        sub = sub.dropna(subset=["numeric_aff"])
        if sub.empty:
            return None, None
        top = sub.sort_values("numeric_aff", ascending=True).iloc[0]
        return top["numeric_aff"], top.get("best_assay_type")

    # ----------------------------------------------------------------
    # Main ligand loop
    # ----------------------------------------------------------------
    for _, comp in df_compound.iterrows():
        ligand_id = comp["molecule_chembl_id"]

        best_name = comp.get("best_name") or comp.get("pref_name")
        smiles = comp.get("canonical_smiles")
        pubchem_cid = comp.get("pubchem_cid")
        cas = comp.get("cas_number")

        # ---------------- Target summary ----------------
        top_target, best_aff, best_assay_type = best_target_for_lig(ligand_id)
        rarg_aff, rarg_assay_type = rarg_affinity_for_lig(ligand_id)
        num_targets = (
            len(targets_by_lig.get_group(ligand_id))
            if ligand_id in targets_by_lig.groups
            else 0
        )

        # ---------------- Pathway summary ----------------
        if ligand_id in paths_by_lig.groups:
            path_df = paths_by_lig.get_group(ligand_id)
            if "pathway_name" in path_df.columns:
                unique_paths = (
                    path_df["pathway_name"]
                    .dropna()
                    .astype(str)
                    .unique()
                    .tolist()
                )
                num_pathways = len(unique_paths)
                pathways_summary = (
                    "; ".join(unique_paths[:3]) if num_pathways > 0 else None
                )
            else:
                num_pathways = len(path_df)
                pathways_summary = None
        else:
            num_pathways = 0
            pathways_summary = None

        # ---------------- Clinical trials ----------------
        if isinstance(clin_by_lig, dict):
            num_clin = 0
            clinical_summary = None
        elif ligand_id in clin_by_lig.groups:
            clin_df = clin_by_lig.get_group(ligand_id)
            num_clin = clin_df["nct_id"].nunique()
            clinical_summary = summarize_clinical_trials(clin_df)
        else:
            num_clin = 0
            clinical_summary = None

        # ---------------- Literature ----------------
        if ligand_id in lit_by_lig.groups:
            lit_df = lit_by_lig.get_group(ligand_id)

            num_pubmed = lit_df[lit_df["source"] == "PubMed"].shape[0]
            num_bindingdb = lit_df[lit_df["source"] == "BindingDB"].shape[0]
            num_chembl = lit_df[lit_df["source"] == "ChEMBL"].shape[0]

            publication_summary = summarize_publications(lit_df)
        else:
            num_pubmed = 0
            num_bindingdb = 0
            num_chembl = 0
            publication_summary = None

        summary_rows.append({
            "ligand_id": ligand_id,
            "best_name": best_name,
            "canonical_smiles": smiles,
            "pubchem_cid": pubchem_cid,
            "cas_number": cas,
            "top_target": top_target,
            "top_target_affinity_nm": best_aff,
            "top_target_assay_type": best_assay_type,
            "rarg_affinity_nm": rarg_aff,
            "rarg_assay_type": rarg_assay_type,
            "num_targets_total": num_targets,
            "num_pathways": num_pathways,
            "num_clinical_trials": num_clin,
            "clinical_trials_summary": clinical_summary,
            "num_pubmed_articles": num_pubmed,
            "num_bindingdb_entries": num_bindingdb,
            "num_chembl_literature": num_chembl,
            "publication_summary": publication_summary,
            "pathways_summary": pathways_summary,
            "notes": "",
        })

    # ----------------------------------------------------------------
    # Save output
    # ----------------------------------------------------------------
    df_summary = pd.DataFrame(summary_rows)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    df_summary.to_csv(output, sep="\t", index=False)

    print(f"\nSummary written: {df_summary.shape} → {output}\n")
    return df_summary


# --------------------------------------------------------------------
# Command-line execution
# --------------------------------------------------------------------
if __name__ == "__main__":
    build_ligand_summary()
