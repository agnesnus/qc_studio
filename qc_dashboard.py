"""
QC Dashboard — Steroid Panel
Levey-Jennings charts with 2SD/3SD bands and outlier flagging.

Run:
    streamlit run qc_dashboard.py --server.address 0.0.0.0
"""

import sqlite3
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path

DB_PATH = Path(__file__).parent / "steroid_panel.db"


@st.cache_data
def get_qc_data():
    conn = sqlite3.connect(str(DB_PATH))
    query = """
        SELECT
            r.run_date,
            a.name AS analyte,
            s.qc_level,
            AVG(res.concentration) AS concentration
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


def flag_outliers(concentrations, sd2_upper, sd2_lower, sd3_upper, sd3_lower):
    flags = [False] * len(concentrations)

    for i, conc in enumerate(concentrations):
        if conc > sd3_upper or conc < sd3_lower:
            flags[i] = True

        if conc > sd2_upper or conc < sd2_lower:
            if i > 0 and (concentrations[i - 1] > sd2_upper or concentrations[i - 1] < sd2_lower):
                flags[i] = True
                flags[i - 1] = True

    return flags


def make_qc_chart(dates, concentrations, mean_val, sd2_upper, sd2_lower, sd3_upper, sd3_lower, title):
    flags = flag_outliers(concentrations, sd2_upper, sd2_lower, sd3_upper, sd3_lower)

    fig = go.Figure()

    # 3SD band
    fig.add_hrect(
        y0=sd3_lower, y1=sd3_upper,
        fillcolor="rgba(255, 99, 71, 0.1)",
        line_width=0,
        annotation_text="±3SD", annotation_position="top left",
    )

    # 2SD band
    fig.add_hrect(
        y0=sd2_lower, y1=sd2_upper,
        fillcolor="rgba(255, 193, 7, 0.15)",
        line_width=0,
        annotation_text="±2SD", annotation_position="top left",
    )

    # Mean line
    fig.add_hline(
        y=mean_val,
        line_dash="dash", line_color="green", line_width=2,
        annotation_text="Mean",
        annotation_position="bottom right",
    )

    # Data line
    fig.add_trace(go.Scatter(
        x=dates, y=concentrations,
        mode="lines+markers",
        marker=dict(size=8, color="royalblue"),
        line=dict(color="royalblue", width=2),
        name="Concentration",
    ))

    # Flagged points
    flagged_dates = [d for d, f in zip(dates, flags) if f]
    flagged_concs = [c for c, f in zip(concentrations, flags) if f]

    if flagged_dates:
        fig.add_trace(go.Scatter(
            x=flagged_dates, y=flagged_concs,
            mode="markers",
            marker=dict(size=12, color="red", symbol="triangle-up", line=dict(width=1, color="darkred")),
            name="Flagged",
        ))

    fig.update_layout(
        title=title,
        xaxis_title="Date",
        yaxis_title="Concentration",
        height=400,
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )

    return fig


def main():
    st.set_page_config(page_title="QC Dashboard — Steroid Panel", layout="wide")
    st.title("QC Dashboard — Steroid Panel")

    if not DB_PATH.exists():
        st.error(f"Database not found at: {DB_PATH}")
        return

    df = get_qc_data()
    if df.empty:
        st.warning("No QC data found in the database.")
        return

    analytes = sorted(df["analyte"].unique())
    selected = st.sidebar.selectbox("Select Hormone", analytes)

    analyte_data = df[df["analyte"] == selected]
    hqc_data = analyte_data[analyte_data["qc_level"] == "High"].reset_index(drop=True)
    lqc_data = analyte_data[analyte_data["qc_level"] == "Low"].reset_index(drop=True)

    for level_name, level_data in [("HQC", hqc_data), ("LQC", lqc_data)]:
        if level_data.empty:
            st.info(f"No {level_name} data for {selected}.")
            continue

        concentrations = level_data["concentration"].tolist()
        dates = [d.replace("-", "/") for d in level_data["run_date"].tolist()]

        mean_val = level_data["concentration"].mean()
        sd = level_data["concentration"].std()

        if pd.isna(sd) or sd == 0:
            st.info(f"Not enough {level_name} data points for {selected} to compute SD.")
            continue

        sd2_upper = mean_val + 2 * sd
        sd2_lower = mean_val - 2 * sd
        sd3_upper = mean_val + 3 * sd
        sd3_lower = mean_val - 3 * sd

        fig = make_qc_chart(
            dates, concentrations,
            mean_val, sd2_upper, sd2_lower, sd3_upper, sd3_lower,
            title=f"{selected} — {level_name}",
        )
        st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
