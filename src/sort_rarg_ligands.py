"""
score_rarg_activity.py

Rank ligands by activity and specificity toward retinoic acid receptor gamma (RARG)
based on a pre-processed activities.tsv file.

"""

import pandas as pd


# -----------------------------------------------------------
# Helper: compute best potency for one ligand + one target
# -----------------------------------------------------------

def _best_pchembl(sub_df, target):
    """Return best pChEMBL value for a target or None."""
    hits = sub_df[sub_df["target_gene_symbol"] == target]
    if hits.empty:
        return None
    return hits["pchembl_value"].max()


# -----------------------------------------------------------
# Helper: compute activity score for RARG potency
# -----------------------------------------------------------

def _score_rarg_activity(best_rarg):
    """Activity score (0–15)."""
    if best_rarg is None:
        return 0
    if best_rarg >= 8.0: return 15
    if best_rarg >= 7.5: return 12
    if best_rarg >= 7.0: return 10
    if best_rarg >= 6.5: return 6
    if best_rarg >= 6.0: return 4
    return 0


# -----------------------------------------------------------
# Helper: compute selectivity score
# -----------------------------------------------------------

def _score_rarg_selectivity(best_rarg, best_rara, best_rarb):
    """Selectivity: reward being stronger on RARG than RARA/RARB."""
    if best_rarg is None:
        return -10  # strong penalty if RARG is missing entirely

    competitors = [x for x in (best_rara, best_rarb) if x is not None]
    if not competitors:
        return 5  # mild reward (no competing data)

    delta = best_rarg - max(competitors)

    if delta >= 1.0: return 20
    if delta >= 0.7: return 12
    if delta >= 0.4: return 6
    if delta >= 0.2: return 3
    if delta >= 0.0: return 0
    return -5  # RARG is weaker than other subtypes


# -----------------------------------------------------------
# Helper: off-target penalty
# -----------------------------------------------------------

def _score_off_target_penalty(sub_df):
    """Penalize the number of non-RAR targets."""
    rar_targets = {"RARG", "RARA", "RARB"}
    all_targets = sub_df["target_gene_symbol"].dropna().unique()
    off_targets = [t for t in all_targets if t not in rar_targets]
    return -2 * len(off_targets), off_targets


# -----------------------------------------------------------
# MAIN FUNCTION
# -----------------------------------------------------------

def score_ligands_by_rarg(activities_tsv, ligands_tsv=None):
    """
    Score ligands for potency and selectivity toward RARG using only
    the three nuclear receptor targets:

        RARA  = CHEMBL2055
        RARB  = CHEMBL2008
        RARG  = CHEMBL2003

    Parameters
    ----------
    activities_tsv : str
        Path to activities.tsv containing columns:
        - molecule_chembl_id
        - target_chembl_id
        - pchembl_value

    ligands_tsv : str or None
        Optional: attach ligand_name or metadata

    Returns
    -------
    pandas.DataFrame
        Columns include:
        - molecule_chembl_id
        - best_RARG
        - best_RARA
        - best_RARB
        - activity_score
        - selectivity_score
        - total_score
    """
    import pandas as pd

    # Load activities file
    acts = pd.read_csv(activities_tsv, sep="\t")
    acts["pchembl_value"] = pd.to_numeric(acts["pchembl_value"], errors="coerce")

    # Required columns check
    required_cols = {"molecule_chembl_id", "target_chembl_id", "pchembl_value"}
    missing = required_cols - set(acts.columns)
    if missing:
        raise ValueError(f"activities.tsv missing required columns: {missing}")

    # Map CHEMBL IDs to canonical names
    TARGET_MAP = {
        "CHEMBL2055": "RARA",
        "CHEMBL2008": "RARB",
        "CHEMBL2003": "RARG",
    }

    # Convert target_chembl_id → gene symbol
    acts["target_gene"] = acts["target_chembl_id"].map(TARGET_MAP)

    def best_p(df, gene):
        sub = df[df["target_gene"] == gene]
        return sub["pchembl_value"].max() if not sub.empty else None

    rows = []
    for chembl_id, sub in acts.groupby("molecule_chembl_id"):

        p_rarg = best_p(sub, "RARG")
        p_rara = best_p(sub, "RARA")
        p_rarb = best_p(sub, "RARB")

        # ----------- RARG Activity -----------
        if p_rarg is None:
            activity_score = 0
        elif p_rarg >= 8.0: activity_score = 20
        elif p_rarg >= 7.5: activity_score = 15
        elif p_rarg >= 7.0: activity_score = 12
        elif p_rarg >= 6.5: activity_score = 8
        elif p_rarg >= 6.0: activity_score = 4
        else: activity_score = 0

        # ----------- Selectivity -------------
        other = [p for p in [p_rara, p_rarb] if p is not None]
        if p_rarg is None:
            selectivity_score = -10
        elif not other:
            selectivity_score = 5
        else:
            delta = p_rarg - max(other)
            if delta >= 1.0: selectivity_score = 20
            elif delta >= 0.7: selectivity_score = 12
            elif delta >= 0.4: selectivity_score = 6
            elif delta >= 0.2: selectivity_score = 3
            elif delta >= 0.0: selectivity_score = 0
            else: selectivity_score = -10

        total_score = activity_score + selectivity_score

        rows.append({
            "molecule_chembl_id": chembl_id,
            "best_RARG": p_rarg,
            "best_RARA": p_rara,
            "best_RARB": p_rarb,
            "activity_score": activity_score,
            "selectivity_score": selectivity_score,
            "total_score": total_score,
        })

    df_out = pd.DataFrame(rows)

    # Optional merge with ligand metadata
    # if ligands_tsv:
    #     ligs = pd.read_csv(ligands_tsv, sep="\t")
    #     if "molecule_chembl_id" in ligs.columns:
    #         df_out = df_out.merge(
    #             ligs[["molecule_chembl_id", "ligand_name"]]
    #             if "ligand_name" in ligs.columns
    #             else ligs[["molecule_chembl_id"]],
    #             on="molecule_chembl_id",
    #             how="left"
    #         )

    df_out = df_out.sort_values("total_score", ascending=False).reset_index(drop=True)
    return df_out




if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Rank ligands by RARG potency + specificity.")
    parser.add_argument("--activities", required=True, help="Path to activities.tsv")
    parser.add_argument("--ligands", required=False, help="Optional path to ligands.tsv")
    parser.add_argument("--id", default="pubchem_cid", help="Identifier column to group by")

    args = parser.parse_args()

    df = score_ligands_by_rarg(args.activities, args.ligands, args.id)
    df.to_csv("ranked_ligands.tsv", sep="\t", index=False)
    print("Wrote ranked_ligands.tsv")
    print(df.head())
