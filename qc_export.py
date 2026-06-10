"""
QC Data Exporter
================
Reads QC_Low and QC_High data from steroid_panel.db and exports
one CSV per hormone with HQC and LQC side by side.

Usage:
    python3 qc_export.py

Outputs:
    Cortisol_QC.csv, Testosterone_QC.csv, etc. on the Desktop.
"""

import sqlite3
import pandas as pd
import numpy as np
from pathlib import Path

DB_PATH = Path(__file__).parent / "steroid_panel.db"
OUTPUT_DIR = Path(__file__).parent


def get_qc_data():
    """Pull all QC results from the database, averaged per run date."""
    conn = sqlite3.connect(str(DB_PATH))
    query = """
        SELECT
            r.run_date,
            a.name as analyte,
            s.qc_level,
            AVG(res.concentration) as concentration
        FROM results res
        JOIN samples s ON res.sample_id = s.sample_id
        JOIN runs r ON s.run_id = r.run_id
        JOIN analytes a ON res.analyte_id = a.analyte_id
        JOIN sample_types st ON s.sample_type_id = st.type_id
        WHERE st.type_code = 'qc'
          AND res.concentration IS NOT NULL
        GROUP BY r.run_date, a.name, s.qc_level
        ORDER BY a.name, r.run_date
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


def format_date(date_str):
    """Convert YYYY-MM-DD to DD/MM/YYYY."""
    parts = date_str.split("-")
    return f"{parts[2]}/{parts[1]}/{parts[0]}"


def export_hormone_csv(analyte_name, hqc_data, lqc_data):
    """Create a CSV for one hormone with HQC and LQC side by side."""
    # Calculate stats for HQC
    hqc_mean = hqc_data["concentration"].mean() if not hqc_data.empty else np.nan
    hqc_sd = hqc_data["concentration"].std() if len(hqc_data) > 1 else np.nan

    # Calculate stats for LQC
    lqc_mean = lqc_data["concentration"].mean() if not lqc_data.empty else np.nan
    lqc_sd = lqc_data["concentration"].std() if len(lqc_data) > 1 else np.nan

    # Get all unique dates across both levels
    all_dates = sorted(set(
        hqc_data["run_date"].tolist() + lqc_data["run_date"].tolist()
    ))

    # Build lookup dicts
    hqc_by_date = dict(zip(hqc_data["run_date"], hqc_data["concentration"]))
    lqc_by_date = dict(zip(lqc_data["run_date"], lqc_data["concentration"]))

    rows = []
    for date in all_dates:
        row = {"Date": format_date(date)}

        # HQC columns
        hqc_conc = hqc_by_date.get(date)
        row["HQC_Conc"] = round(hqc_conc, 4) if hqc_conc is not None else ""
        row["HQC_Mean"] = round(hqc_mean, 4) if not np.isnan(hqc_mean) else ""
        row["HQC_+2SD"] = round(hqc_mean + 2 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""
        row["HQC_-2SD"] = round(hqc_mean - 2 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""
        row["HQC_+3SD"] = round(hqc_mean + 3 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""
        row["HQC_-3SD"] = round(hqc_mean - 3 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""

        # LQC columns
        lqc_conc = lqc_by_date.get(date)
        row["LQC_Conc"] = round(lqc_conc, 4) if lqc_conc is not None else ""
        row["LQC_Mean"] = round(lqc_mean, 4) if not np.isnan(lqc_mean) else ""
        row["LQC_+2SD"] = round(lqc_mean + 2 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""
        row["LQC_-2SD"] = round(lqc_mean - 2 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""
        row["LQC_+3SD"] = round(lqc_mean + 3 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""
        row["LQC_-3SD"] = round(lqc_mean - 3 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""

        rows.append(row)

    df_out = pd.DataFrame(rows)
    filename = f"{analyte_name}_QC.csv"
    filepath = OUTPUT_DIR / filename
    df_out.to_csv(filepath, index=False)
    return filepath


def main():
    if not DB_PATH.exists():
        print(f"Database not found: {DB_PATH}")
        print("Run 'python3 steroid_panel_lims.py init' and import data first.")
        return

    df = get_qc_data()
    if df.empty:
        print("No QC data found in the database.")
        return

    analytes = sorted(df["analyte"].unique())
    exported = []

    for analyte in analytes:
        analyte_data = df[df["analyte"] == analyte]
        hqc_data = analyte_data[analyte_data["qc_level"] == "High"].reset_index(drop=True)
        lqc_data = analyte_data[analyte_data["qc_level"] == "Low"].reset_index(drop=True)

        if hqc_data.empty and lqc_data.empty:
            continue

        filepath = export_hormone_csv(analyte, hqc_data, lqc_data)
        exported.append(filepath.name)

    print(f"Exported {len(exported)} CSV files to {OUTPUT_DIR}:")
    for name in exported:
        print(f"  {name}")


if __name__ == "__main__":
    main()
