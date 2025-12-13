import pandas as pd

def normalize_phase(val):
    """
    Convert Open Targets 'phase' strings into sortable numeric values.
    Examples:
        'Phase 4' -> 4
        'Phase 3' -> 3
        'Phase 2/3' -> 2.5
        'Preclinical' -> 0
        None / '' -> -1
    """
    if val is None or pd.isna(val):
        return -1

    text = str(val).strip().lower()

    if text.startswith("phase"):
        # handle things like "Phase 2/3"
        parts = text.replace("phase", "").strip()
        if "/" in parts:
            a, b = parts.split("/", 1)
            try:
                return (float(a) + float(b)) / 2.0
            except:
                return -1
        try:
            return float(parts)
        except:
            return -1

    if "preclinical" in text:
        return 0

    return -1  # fallback


def drugs_with_highest_phase(df):
    """
    Given df_ot_drugs, return each unique drug with its highest phase.
    """
    # Normalize clinical phase into numeric
    df = df.copy()
    df["phase_numeric"] = df["phase"].apply(normalize_phase)

    # Pick the row with the highest phase per drug
    idx = df.groupby("drug_name")["phase_numeric"].idxmax()
    top_phase_rows = df.loc[idx].sort_values("phase_numeric", ascending=False)

    return top_phase_rows[[
        "drug_name",
        "drug_id",
        "phase",
        "phase_numeric",
        "status",
        "indication_name",
        "is_approved"
    ]]