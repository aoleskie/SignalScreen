"""
sort_rarg_ligands.py

Rank ligands by activity and specificity toward retinoic acid receptor gamma (RARG)
based on a pre-processed activities.tsv file.

Scoring is continuous rather than binned. The previous binned version saturated at
20 + 20 = 40, which produced a 10-way tie at the top of the table and made the
"top N" cut an arbitrary slice of that tie.
"""

import pandas as pd

# RAR subtype ChEMBL target IDs
TARGET_MAP = {
    "CHEMBL2055": "RARA",
    "CHEMBL2008": "RARB",
    "CHEMBL2003": "RARG",
}

# Activity: linear in pChEMBL, 0 at ACTIVITY_FLOOR, 20 at ACTIVITY_CEIL.
# The floor keeps the old "below pChEMBL 6 is worthless" cutoff; the ceiling sits
# above the strongest value in the dataset so nothing saturates.
ACTIVITY_FLOOR = 6.0
ACTIVITY_CEIL = 10.5
ACTIVITY_MAX = 20.0

# Selectivity: linear in delta = pRARG - max(pRARA, pRARB).
SELECTIVITY_PER_LOG = 10.0
SELECTIVITY_MAX = 20.0
SELECTIVITY_MIN = -10.0

# Applied when RARG potency is known but neither other subtype was measured.
# Unmeasured is not the same as selective, so this is a mild reward only; the
# selectivity_evidence column records which ligands got it.
SELECTIVITY_NO_DATA = 5.0


def _score_rarg_activity(p_rarg):
    """Potency score, 0-20, linear in pChEMBL between the floor and ceiling."""
    if p_rarg is None:
        return 0.0
    scaled = ACTIVITY_MAX * (p_rarg - ACTIVITY_FLOOR) / (ACTIVITY_CEIL - ACTIVITY_FLOOR)
    return round(min(max(scaled, 0.0), ACTIVITY_MAX), 3)


def _score_rarg_selectivity(p_rarg, p_rara, p_rarb):
    """Selectivity score, -10 to 20, linear in the log-unit gap over RARA/RARB."""
    if p_rarg is None:
        return SELECTIVITY_MIN, None, "no_rarg_data"

    competitors = [p for p in (p_rara, p_rarb) if p is not None]
    if not competitors:
        return SELECTIVITY_NO_DATA, None, "no_subtype_data"

    delta = p_rarg - max(competitors)
    scaled = SELECTIVITY_PER_LOG * delta
    score = min(max(scaled, SELECTIVITY_MIN), SELECTIVITY_MAX)
    return round(score, 3), round(delta, 3), "measured"


def score_ligands_by_rarg(activities_tsv, ligands_tsv=None):
    """
    Score ligands for potency and selectivity toward RARG using direct binding
    data against the three RAR subtypes (RARA/RARB/RARG).

    Parameters
    ----------
    activities_tsv : str or Path
        Path to activities.tsv with columns molecule_chembl_id, target_chembl_id,
        pchembl_value.
    ligands_tsv : str or Path or None
        Optional ligands.tsv; when given, canonical_smiles and names are merged in.

    Returns
    -------
    pandas.DataFrame
        Sorted best-first on (total_score, best_RARG, selectivity_delta).
    """
    acts = pd.read_csv(activities_tsv, sep="\t")

    required_cols = {"molecule_chembl_id", "target_chembl_id", "pchembl_value"}
    missing = required_cols - set(acts.columns)
    if missing:
        raise ValueError(f"activities.tsv missing required columns: {missing}")

    acts["pchembl_value"] = pd.to_numeric(acts["pchembl_value"], errors="coerce")
    acts["target_gene"] = acts["target_chembl_id"].map(TARGET_MAP)

    def best_p(df, gene):
        sub = df[df["target_gene"] == gene]
        val = sub["pchembl_value"].max()
        return None if pd.isna(val) else float(val)

    rows = []
    for chembl_id, sub in acts.groupby("molecule_chembl_id"):
        p_rarg = best_p(sub, "RARG")
        p_rara = best_p(sub, "RARA")
        p_rarb = best_p(sub, "RARB")

        activity_score = _score_rarg_activity(p_rarg)
        selectivity_score, delta, evidence = _score_rarg_selectivity(p_rarg, p_rara, p_rarb)

        rows.append({
            "molecule_chembl_id": chembl_id,
            "best_RARG": p_rarg,
            "best_RARA": p_rara,
            "best_RARB": p_rarb,
            "selectivity_delta": delta,
            "selectivity_evidence": evidence,
            "activity_score": activity_score,
            "selectivity_score": selectivity_score,
            "total_score": round(activity_score + selectivity_score, 3),
        })

    df_out = pd.DataFrame(rows)

    # Explicit tie-break: raw RARG potency first, then the selectivity gap.
    df_out = df_out.sort_values(
        ["total_score", "best_RARG", "selectivity_delta"],
        ascending=False,
        na_position="last",
    ).reset_index(drop=True)

    if ligands_tsv:
        ligs = pd.read_csv(ligands_tsv, sep="\t")
        keep = [c for c in ("molecule_chembl_id", "canonical_smiles", "names") if c in ligs.columns]
        if len(keep) > 1:
            df_out = df_out.merge(ligs[keep], on="molecule_chembl_id", how="left")

    return df_out


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Rank ligands by RARG potency + specificity.")
    parser.add_argument("--activities", required=True, help="Path to activities.tsv")
    parser.add_argument("--ligands", required=False, help="Optional path to ligands.tsv")
    parser.add_argument("--out", default="data/processed/sorted_ligands_RARG.tsv",
                        help="Output TSV path")
    args = parser.parse_args()

    df = score_ligands_by_rarg(args.activities, args.ligands)
    from pathlib import Path
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {len(df)} scored ligands to {args.out}")
    print(df.head().to_string())
