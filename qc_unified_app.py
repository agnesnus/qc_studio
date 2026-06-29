"""
QC Studio — Unified Application
================================
Integrated platform for steroid panel database management, QC data export, and dashboard visualization.

Features:
1. Steroid Panel Database: SQLite database for LC-MS/MS steroid panel results (built from uploaded files)
2. QC Export: Export CSV files with HQC and LQC values for all hormones
3. QC Dashboard: Interactive Levey-Jennings charts with 2SD/3SD bands

Run:
    streamlit run qc_unified_app.py --server.address localhost --server.port 8501
"""

import sqlite3
import re
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from pathlib import Path
from dataclasses import dataclass
from typing import Optional
from datetime import datetime
import tempfile
import os
from io import BytesIO
from urllib.parse import quote

# ==============================================================================
# CONFIGURATION
# ==============================================================================

DB_PATH = Path(__file__).parent / "steroid_panel.db"

ANALYTES = [
    {"name": "11-deoxycorticosterone", "panel": 1, "display_order": 1},
    {"name": "11-deoxycortisol", "panel": 1, "display_order": 2},
    {"name": "17-hydroxyprogesterone", "panel": 1, "display_order": 3},
    {"name": "21-deoxycortisol", "panel": 1, "display_order": 4},
    {"name": "Aldosterone", "panel": 1, "display_order": 5},
    {"name": "Androstenedione", "panel": 1, "display_order": 6},
    {"name": "Corticosterone", "panel": 1, "display_order": 7},
    {"name": "Cortisol", "panel": 1, "display_order": 8},
    {"name": "Cortisone", "panel": 1, "display_order": 9},
    {"name": "Dexamethasone", "panel": 1, "display_order": 10},
    {"name": "DHEA", "panel": 1, "display_order": 11},
    {"name": "DHEAS", "panel": 1, "display_order": 12},
    {"name": "Dihydrotestosterone", "panel": 1, "display_order": 13},
    {"name": "Progesterone", "panel": 1, "display_order": 14},
    {"name": "Testosterone", "panel": 1, "display_order": 15},
    {"name": "17-hydroxypregnenolone", "panel": 2, "display_order": 16},
    {"name": "Estradiol", "panel": 2, "display_order": 17},
    {"name": "Estrone", "panel": 2, "display_order": 18},
]

SAMPLE_TYPES = [
    {"type_code": "calibrator", "description": "Calibration standards (Cal 0 through Cal F)"},
    {"type_code": "qc", "description": "Quality control samples (Low/High)"},
    {"type_code": "patient", "description": "Patient specimens"},
    {"type_code": "eqa", "description": "External quality assessment / proficiency testing"},
    {"type_code": "blank", "description": "Solvent blanks"},
    {"type_code": "process_blank", "description": "Process/extraction blanks"},
]

# ==============================================================================
# SQL SCHEMA
# ==============================================================================

SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS runs (
    run_id          INTEGER PRIMARY KEY AUTOINCREMENT,
    run_date        TEXT NOT NULL,
    panel           INTEGER NOT NULL,
    source_filename TEXT NOT NULL UNIQUE,
    method_name     TEXT,
    data_path       TEXT,
    uploaded_by     TEXT,
    imported_at     TEXT NOT NULL DEFAULT (datetime('now')),
    notes           TEXT
);

CREATE TABLE IF NOT EXISTS analytes (
    analyte_id      INTEGER PRIMARY KEY AUTOINCREMENT,
    name            TEXT NOT NULL UNIQUE,
    panel           INTEGER NOT NULL,
    display_order   INTEGER
);

CREATE TABLE IF NOT EXISTS sample_types (
    type_id         INTEGER PRIMARY KEY AUTOINCREMENT,
    type_code       TEXT NOT NULL UNIQUE,
    description     TEXT
);

CREATE TABLE IF NOT EXISTS samples (
    sample_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id          INTEGER NOT NULL REFERENCES runs(run_id),
    data_filename   TEXT NOT NULL,
    sample_name     TEXT,
    sample_type_id  INTEGER NOT NULL REFERENCES sample_types(type_id),
    instrument_type TEXT,
    acquisition_datetime TEXT,
    autosampler_position TEXT,
    sample_group    TEXT,
    collection_date TEXT,
    patient_sequence TEXT,
    calibrator_level TEXT,
    qc_level        TEXT,
    qc_replicate    INTEGER,
    eqa_scheme      TEXT,
    eqa_year        INTEGER,
    eqa_round       INTEGER,
    eqa_sample_code TEXT,
    eqa_replicate   TEXT,
    UNIQUE(run_id, data_filename)
);

CREATE TABLE IF NOT EXISTS results (
    result_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id       INTEGER NOT NULL REFERENCES samples(sample_id),
    analyte_id      INTEGER NOT NULL REFERENCES analytes(analyte_id),
    concentration   REAL,
    UNIQUE(sample_id, analyte_id)
);

CREATE TABLE IF NOT EXISTS qc_targets (
    target_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    analyte_id      INTEGER NOT NULL REFERENCES analytes(analyte_id),
    qc_level        TEXT NOT NULL,
    lot_number      TEXT,
    target_mean     REAL NOT NULL,
    target_sd       REAL NOT NULL,
    effective_from  TEXT NOT NULL,
    effective_to    TEXT,
    UNIQUE(analyte_id, qc_level, lot_number, effective_from)
);

CREATE TABLE IF NOT EXISTS eqa_targets (
    target_id       INTEGER PRIMARY KEY AUTOINCREMENT,
    analyte_id      INTEGER NOT NULL REFERENCES analytes(analyte_id),
    scheme          TEXT NOT NULL,
    year            INTEGER NOT NULL,
    round           INTEGER NOT NULL,
    sample_code     TEXT NOT NULL,
    consensus_mean  REAL,
    consensus_sd    REAL,
    UNIQUE(analyte_id, scheme, year, round, sample_code)
);

CREATE INDEX IF NOT EXISTS idx_results_analyte ON results(analyte_id);
CREATE INDEX IF NOT EXISTS idx_results_sample ON results(sample_id);
CREATE INDEX IF NOT EXISTS idx_samples_type ON samples(sample_type_id);
CREATE INDEX IF NOT EXISTS idx_samples_run ON samples(run_id);
CREATE INDEX IF NOT EXISTS idx_runs_date ON runs(run_date);
CREATE INDEX IF NOT EXISTS idx_samples_qc ON samples(qc_level) WHERE qc_level IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_samples_eqa ON samples(eqa_scheme, eqa_year, eqa_round, eqa_sample_code)
    WHERE eqa_scheme IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_qc_targets_lookup ON qc_targets(analyte_id, qc_level, effective_from);
"""

# ==============================================================================
# SAMPLE CLASSIFIER
# ==============================================================================

@dataclass
class SampleInfo:
    data_filename: str
    sample_type: str
    calibrator_level: Optional[str] = None
    qc_level: Optional[str] = None
    qc_replicate: Optional[int] = None
    collection_date: Optional[str] = None
    patient_sequence: Optional[str] = None
    eqa_scheme: Optional[str] = None
    eqa_year: Optional[int] = None
    eqa_round: Optional[int] = None
    eqa_sample_code: Optional[str] = None
    eqa_replicate: Optional[str] = None


def classify_sample(data_filename: str) -> SampleInfo:
    """Parse data filename into structured sample information."""
    base = re.sub(r"(_P[12])?\.d$", "", data_filename)

    cal_match = re.match(r"^Cal\s+([0A-F])$", base)
    if cal_match:
        return SampleInfo(data_filename=data_filename, sample_type="calibrator", calibrator_level=cal_match.group(1))

    qc_match = re.match(r"^QC_(Low|High)(\d+)$", base)
    if qc_match:
        return SampleInfo(
            data_filename=data_filename, sample_type="qc",
            qc_level=qc_match.group(1), qc_replicate=int(qc_match.group(2))
        )

    if re.match(r"^Blank\d*$", base):
        return SampleInfo(data_filename=data_filename, sample_type="blank")

    if re.match(r"^(PBlank|PB)\d*$", base):
        return SampleInfo(data_filename=data_filename, sample_type="process_blank")

    eqa_match = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])-([a-z])$", base)
    if eqa_match:
        return SampleInfo(
            data_filename=data_filename, sample_type="eqa",
            eqa_scheme=eqa_match.group(1), eqa_year=int(eqa_match.group(2)),
            eqa_round=int(eqa_match.group(3)), eqa_sample_code=eqa_match.group(4),
            eqa_replicate=eqa_match.group(5)
        )

    eqa_special = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])_(.+)$", base)
    if eqa_special:
        return SampleInfo(
            data_filename=data_filename, sample_type="eqa",
            eqa_scheme=eqa_special.group(1), eqa_year=int(eqa_special.group(2)),
            eqa_round=int(eqa_special.group(3)), eqa_sample_code=eqa_special.group(4),
            eqa_replicate=eqa_special.group(5)
        )

    eqa_new_rep = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])([a-z])$", base)
    if eqa_new_rep:
        return SampleInfo(
            data_filename=data_filename, sample_type="eqa",
            eqa_scheme=eqa_new_rep.group(1), eqa_year=int(eqa_new_rep.group(2)),
            eqa_round=int(eqa_new_rep.group(3)), eqa_sample_code=eqa_new_rep.group(4),
            eqa_replicate=eqa_new_rep.group(5)
        )

    eqa_new = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])$", base)
    if eqa_new:
        return SampleInfo(
            data_filename=data_filename, sample_type="eqa",
            eqa_scheme=eqa_new.group(1), eqa_year=int(eqa_new.group(2)),
            eqa_round=int(eqa_new.group(3)), eqa_sample_code=eqa_new.group(4), eqa_replicate=None
        )

    patient_match = re.match(r"^(\d{8})_(\w+)$", base)
    if patient_match:
        date_str = patient_match.group(1)
        formatted_date = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
        return SampleInfo(
            data_filename=data_filename, sample_type="patient",
            collection_date=formatted_date, patient_sequence=patient_match.group(2)
        )

    return SampleInfo(data_filename=data_filename, sample_type="patient")


def classify_from_instrument_type(instrument_type: str, level: str, data_filename: str) -> SampleInfo:
    """Use instrument-assigned Type column to classify."""
    itype = str(instrument_type).strip() if pd.notna(instrument_type) else ""
    lvl = str(level).strip() if pd.notna(level) else ""

    if itype == "DoubleBlank":
        return SampleInfo(data_filename=data_filename, sample_type="blank")

    if itype == "Blank":
        return SampleInfo(data_filename=data_filename, sample_type="process_blank")

    if itype == "MatrixBlank":
        return SampleInfo(data_filename=data_filename, sample_type="calibrator", calibrator_level="0")

    if itype == "Cal":
        return SampleInfo(data_filename=data_filename, sample_type="calibrator", calibrator_level=lvl)

    if itype == "QC":
        qc_match = re.match(r"^QC_(Low|High)(\d+)", re.sub(r"(_P[12])?\.d$", "", data_filename))
        qc_level = lvl if lvl else None
        qc_replicate = int(qc_match.group(2)) if qc_match else 1
        if qc_match and not qc_level:
            qc_level = qc_match.group(1)
        return SampleInfo(
            data_filename=data_filename, sample_type="qc",
            qc_level=qc_level, qc_replicate=qc_replicate
        )

    return classify_sample(data_filename)


# ==============================================================================
# DATABASE OPERATIONS
# ==============================================================================

def ensure_uploaded_by_column(conn):
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='runs'")
    if cursor.fetchone():
        cursor.execute("PRAGMA table_info(runs)")
        columns = [row[1] for row in cursor.fetchall()]
        if "uploaded_by" not in columns:
            cursor.execute("ALTER TABLE runs ADD COLUMN uploaded_by TEXT")


def get_connection(db_path=None):
    db_path = db_path or DB_PATH
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    ensure_uploaded_by_column(conn)
    return conn


def ensure_db_initialized(db_path=None):
    """Create schema and seed reference data on-demand (from first file upload)."""
    conn = get_connection(db_path)
    cursor = conn.cursor()
    
    # Create schema if not exists
    conn.executescript(SCHEMA_SQL)

    # Seed analytes if not already present
    cursor.execute("SELECT COUNT(*) FROM analytes")
    if cursor.fetchone()[0] == 0:
        for analyte in ANALYTES:
            conn.execute(
                "INSERT OR IGNORE INTO analytes (name, panel, display_order) VALUES (?, ?, ?)",
                (analyte["name"], analyte["panel"], analyte["display_order"])
            )

    # Seed sample types if not already present
    cursor.execute("SELECT COUNT(*) FROM sample_types")
    if cursor.fetchone()[0] == 0:
        for st_item in SAMPLE_TYPES:
            conn.execute(
                "INSERT OR IGNORE INTO sample_types (type_code, description) VALUES (?, ?)",
                (st_item["type_code"], st_item["description"])
            )

    conn.commit()
    conn.close()


# ==============================================================================
# FORMAT DETECTION AND FILENAME PARSING
# ==============================================================================

def detect_format(csv_path: str) -> str:
    """Detect whether CSV is old or new format by checking header row."""
    df_header = pd.read_csv(csv_path, nrows=1, header=0)
    second_row = df_header.iloc[0].tolist() if len(df_header) > 0 else []
    second_row_strs = [str(v).strip() for v in second_row if pd.notna(v)]
    if "Data Path" in second_row_strs:
        return "new"
    return "old"


def parse_filename_old(filepath: str) -> dict:
    """Parse '20260206_Panel1_conc(in).csv' → run metadata."""
    basename = Path(filepath).name
    match = re.match(r"^(\d{8})_Panel(\d+)_conc\(in\)\.csv$", basename)
    if not match:
        raise ValueError(f"Filename does not match old pattern: {basename}")
    date_str = match.group(1)
    run_date = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
    panel = int(match.group(2))
    return {"run_date": run_date, "panel": panel, "source_filename": basename, "method_name": None}


def parse_filename_new(filepath: str) -> dict:
    """Parse 'ISOSP-23_20260420_P1(Sheet1).csv' → run metadata."""
    basename = Path(filepath).name
    match = re.match(r"^(.+?)_(\d{8})_P(\d+)(?:-\w+)?\(Sheet\d+\)\.csv$", basename)
    if not match:
        raise ValueError(f"Filename does not match new pattern: {basename}")
    method_name = match.group(1)
    date_str = match.group(2)
    run_date = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
    panel = int(match.group(3))
    return {"run_date": run_date, "panel": panel, "source_filename": basename, "method_name": method_name}


# ==============================================================================
# CSV IMPORTER — OLD FORMAT
# ==============================================================================

def import_csv_old(csv_path: str, db_path=None, uploaded_by=None):
    """Import old-format CSV."""
    ensure_db_initialized(db_path)
    conn = get_connection(db_path)
    cursor = conn.cursor()

    uploaded_by = str(uploaded_by).strip().upper() if uploaded_by else None
    meta = parse_filename_old(csv_path)

    cursor.execute("SELECT run_id FROM runs WHERE source_filename = ?", (meta["source_filename"],))
    if cursor.fetchone():
        return f"Already imported: {meta['source_filename']}"

    cursor.execute(
        "INSERT INTO runs (run_date, panel, source_filename, method_name, uploaded_by) VALUES (?, ?, ?, ?, ?)",
        (meta["run_date"], meta["panel"], meta["source_filename"], meta["method_name"], uploaded_by)
    )
    run_id = cursor.lastrowid

    df = pd.read_csv(csv_path, header=0)
    df = df.iloc[1:].reset_index(drop=True)

    analyte_columns = [col.replace(" Results", "").strip() for col in df.columns[1:]]

    analyte_id_map = {}
    for name in analyte_columns:
        cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (name,))
        row = cursor.fetchone()
        if row:
            analyte_id_map[name] = row[0]

    cursor.execute("SELECT type_code, type_id FROM sample_types")
    type_map = dict(cursor.fetchall())

    imported_count = 0
    for _, row in df.iterrows():
        data_filename = str(row.iloc[0]).strip()
        if not data_filename or data_filename == "nan":
            continue

        info = classify_sample(data_filename)
        sample_type_id = type_map[info.sample_type]

        cursor.execute(
            """INSERT INTO samples (
                run_id, data_filename, sample_name, sample_type_id, instrument_type,
                acquisition_datetime, autosampler_position, sample_group,
                collection_date, patient_sequence,
                calibrator_level, qc_level, qc_replicate,
                eqa_scheme, eqa_year, eqa_round, eqa_sample_code, eqa_replicate
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (run_id, data_filename, None, sample_type_id, None, None, None, None,
             info.collection_date, info.patient_sequence, info.calibrator_level, info.qc_level, info.qc_replicate,
             info.eqa_scheme, info.eqa_year, info.eqa_round, info.eqa_sample_code, info.eqa_replicate)
        )
        sample_id = cursor.lastrowid

        for i, analyte_name in enumerate(analyte_columns):
            raw_value = row.iloc[i + 1]
            concentration = None
            if pd.notna(raw_value) and str(raw_value).strip() != "":
                try:
                    concentration = float(raw_value)
                except ValueError:
                    pass

            cursor.execute(
                "INSERT INTO results (sample_id, analyte_id, concentration) VALUES (?, ?, ?)",
                (sample_id, analyte_id_map[analyte_name], concentration)
            )

        imported_count += 1

    conn.commit()
    conn.close()
    return f"Imported {imported_count} samples from {meta['source_filename']}"


# ==============================================================================
# CSV IMPORTER — NEW FORMAT
# ==============================================================================

def import_csv_new(csv_path: str, db_path=None, uploaded_by=None):
    """Import new-format CSV with full metadata columns."""
    ensure_db_initialized(db_path)
    conn = get_connection(db_path)
    cursor = conn.cursor()

    uploaded_by = str(uploaded_by).strip().upper() if uploaded_by else None
    meta = parse_filename_new(csv_path)

    cursor.execute("SELECT run_id FROM runs WHERE source_filename = ?", (meta["source_filename"],))
    if cursor.fetchone():
        return f"Already imported: {meta['source_filename']}"

    df_raw = pd.read_csv(csv_path, header=None)

    top_headers = df_raw.iloc[0].tolist()

    analyte_start_idx = None
    analyte_columns = []
    for i, h in enumerate(top_headers):
        h_str = str(h).strip() if pd.notna(h) else ""
        if "Results" in h_str:
            if analyte_start_idx is None:
                analyte_start_idx = i
            analyte_columns.append(h_str.replace(" Results", "").strip())

    if analyte_start_idx is None:
        return "Error: Could not find analyte Results columns in CSV header"

    sub_headers = df_raw.iloc[1].tolist()
    meta_col_names = [str(v).strip() if pd.notna(v) else "" for v in sub_headers[:analyte_start_idx]]

    def find_col(name):
        for i, c in enumerate(meta_col_names):
            if c == name:
                return i
        return None

    col_name = find_col("Name")
    col_data_file = find_col("Data File")
    col_data_path = find_col("Data Path")
    col_type = find_col("Type")
    col_level = find_col("Level")
    col_acq_datetime = find_col("Acq. Date-Time")
    col_sample_group = find_col("Sample Group")
    col_pos = find_col("Pos.")

    first_data_path = None
    if col_data_path is not None and len(df_raw) > 2:
        first_data_path = str(df_raw.iloc[2, col_data_path]).strip() if pd.notna(df_raw.iloc[2, col_data_path]) else None

    cursor.execute(
        "INSERT INTO runs (run_date, panel, source_filename, method_name, data_path, uploaded_by) VALUES (?, ?, ?, ?, ?, ?)",
        (meta["run_date"], meta["panel"], meta["source_filename"], meta["method_name"], first_data_path, uploaded_by)
    )
    run_id = cursor.lastrowid

    analyte_id_map = {}
    for name in analyte_columns:
        cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (name,))
        row = cursor.fetchone()
        if row:
            analyte_id_map[name] = row[0]

    cursor.execute("SELECT type_code, type_id FROM sample_types")
    type_map = dict(cursor.fetchall())

    imported_count = 0
    for row_idx in range(2, len(df_raw)):
        row = df_raw.iloc[row_idx]

        data_filename = None
        if col_data_file is not None:
            val = row.iloc[col_data_file]
            data_filename = str(val).strip() if pd.notna(val) else None

        if not data_filename or data_filename == "nan":
            continue

        sample_name = str(row.iloc[col_name]).strip() if col_name is not None and pd.notna(row.iloc[col_name]) else None
        instrument_type = str(row.iloc[col_type]).strip() if col_type is not None and pd.notna(row.iloc[col_type]) else None
        level = str(row.iloc[col_level]).strip() if col_level is not None and pd.notna(row.iloc[col_level]) else None
        acq_datetime = str(row.iloc[col_acq_datetime]).strip() if col_acq_datetime is not None and pd.notna(row.iloc[col_acq_datetime]) else None
        sample_group = str(row.iloc[col_sample_group]).strip() if col_sample_group is not None and pd.notna(row.iloc[col_sample_group]) else None
        position = str(row.iloc[col_pos]).strip() if col_pos is not None and pd.notna(row.iloc[col_pos]) else None

        if instrument_type and instrument_type != "nan":
            info = classify_from_instrument_type(instrument_type, level, data_filename)
        else:
            info = classify_sample(data_filename)

        sample_type_id = type_map[info.sample_type]

        cursor.execute(
            """INSERT INTO samples (
                run_id, data_filename, sample_name, sample_type_id, instrument_type,
                acquisition_datetime, autosampler_position, sample_group,
                collection_date, patient_sequence,
                calibrator_level, qc_level, qc_replicate,
                eqa_scheme, eqa_year, eqa_round, eqa_sample_code, eqa_replicate
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (run_id, data_filename, sample_name, sample_type_id, instrument_type,
             acq_datetime, position, sample_group, info.collection_date, info.patient_sequence,
             info.calibrator_level, info.qc_level, info.qc_replicate,
             info.eqa_scheme, info.eqa_year, info.eqa_round, info.eqa_sample_code, info.eqa_replicate)
        )
        sample_id = cursor.lastrowid

        for i, analyte_name in enumerate(analyte_columns):
            col_idx = analyte_start_idx + i
            raw_value = row.iloc[col_idx] if col_idx < len(row) else None
            concentration = None
            if pd.notna(raw_value) and str(raw_value).strip() != "":
                try:
                    concentration = float(raw_value)
                except ValueError:
                    pass

            cursor.execute(
                "INSERT INTO results (sample_id, analyte_id, concentration) VALUES (?, ?, ?)",
                (sample_id, analyte_id_map[analyte_name], concentration)
            )

        imported_count += 1

    conn.commit()
    conn.close()
    return f"Imported {imported_count} samples from {meta['source_filename']}"


def import_csv(csv_path: str, db_path=None, uploaded_by=None):
    """Auto-detect format and import CSV."""
    fmt = detect_format(csv_path)
    if fmt == "old":
        return import_csv_old(csv_path, db_path, uploaded_by)
    else:
        return import_csv_new(csv_path, db_path, uploaded_by)


def normalize_qc_level(value):
    if pd.isna(value):
        return None

    label = str(value).strip().lower()
    if not label:
        return None

    if label in {"hqc", "high", "high qc", "high control", "high_level", "high level", "highqc"}:
        return "High"
    if label in {"lqc", "low", "low qc", "low control", "low_level", "low level", "lowqc"}:
        return "Low"

    if "high" in label and "low" not in label:
        return "High"
    if "low" in label and "high" not in label:
        return "Low"

    if label == "h":
        return "High"
    if label == "l":
        return "Low"

    return None


def find_column(df, tokens):
    for col in df.columns:
        lower = str(col).strip().lower()
        for token in tokens:
            if token in lower:
                return col
    return None


def parse_date_value(value):
    if pd.isna(value):
        return None

    try:
        dt = pd.to_datetime(str(value), dayfirst=False, errors='coerce')
    except Exception:
        return None

    if pd.isna(dt):
        return None
    return dt.strftime("%Y-%m-%d")


def extract_date_from_filename(filename):
    if not filename:
        return None

    text = str(filename)
    match = re.search(r"(\d{4}-\d{2}-\d{2})", text)
    if match:
        return match.group(1)

    match = re.search(r"(\d{8})", text)
    if match:
        date_str = match.group(1)
        return f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"

    return None


def normalize_analyte_name(name):
    if pd.isna(name):
        return ""
    normalized = re.sub(r"[^a-z0-9]", "", str(name).strip().lower())
    return normalized


def find_analyte_name_in_workbook(sheet_name, sample_row=None):
    analyte_map = {normalize_analyte_name(a["name"]): a["name"] for a in ANALYTES}
    normalized_sheet = normalize_analyte_name(sheet_name)
    if normalized_sheet in analyte_map:
        return analyte_map[normalized_sheet]

    if sample_row is not None:
        for value in sample_row:
            if pd.notna(value):
                normalized_value = normalize_analyte_name(value)
                if normalized_value in analyte_map:
                    return analyte_map[normalized_value]

    return None


def find_qc_summary_header_row(df):
    for row_idx in range(min(len(df) - 1, 40)):
        row_vec = df.iloc[row_idx].tolist()
        row_strs = [str(v).strip().lower() if pd.notna(v) else "" for v in row_vec]
        if any("hqc" in s for s in row_strs) and any("lqc" in s for s in row_strs) and any(
            token in s for s in row_strs for token in ["qc mean", "%cv", "+2sd", "-2sd"]
        ):
            return row_idx, row_strs
    return None, None


def find_header_index(row_strs, tokens, start=0, end=None):
    end = end if end is not None else len(row_strs)
    for idx in range(start, end):
        cell = row_strs[idx]
        if not cell:
            continue
        if all(token in cell for token in tokens):
            return idx
    return None


def calculate_sd(mean, upper2=None, lower2=None, upper3=None, lower3=None, cv=None):
    if mean is None or pd.isna(mean):
        return None

    if upper2 is not None and not pd.isna(upper2):
        try:
            return abs(float(upper2) - float(mean)) / 2
        except Exception:
            pass

    if lower2 is not None and not pd.isna(lower2):
        try:
            return abs(float(mean) - float(lower2)) / 2
        except Exception:
            pass

    if upper3 is not None and not pd.isna(upper3):
        try:
            return abs(float(upper3) - float(mean)) / 3
        except Exception:
            pass

    if lower3 is not None and not pd.isna(lower3):
        try:
            return abs(float(mean) - float(lower3)) / 3
        except Exception:
            pass

    if cv is not None and not pd.isna(cv):
        try:
            return float(mean) * float(cv) / 100.0
        except Exception:
            pass

    return None


def parse_qc_targets_from_sheet(sheet_name, df, file_date):
    analyte = find_analyte_name_in_workbook(sheet_name, sample_row=df.iloc[0].tolist() if len(df) > 0 else None)
    if analyte is None:
        return []

    header_row_idx, header_row = find_qc_summary_header_row(df)
    if header_row_idx is None or header_row_idx + 1 >= len(df):
        return []

    value_row = df.iloc[header_row_idx + 1].tolist()
    hqc_start = next((i for i, cell in enumerate(header_row) if "hqc" in cell and "%cv" in cell), None)
    lqc_start = next((i for i, cell in enumerate(header_row) if "lqc" in cell and "%cv" in cell), None)
    if hqc_start is None or lqc_start is None:
        return []

    hqc_end = lqc_start
    lqc_end = len(header_row)

    hqc_mean_col = find_header_index(header_row, ["qc mean"], start=hqc_start, end=hqc_end)
    hqc_upper2_col = find_header_index(header_row, ["+2sd"], start=hqc_start, end=hqc_end)
    hqc_lower2_col = find_header_index(header_row, ["-2sd"], start=hqc_start, end=hqc_end)
    hqc_upper3_col = find_header_index(header_row, ["+3sd"], start=hqc_start, end=hqc_end)
    hqc_lower3_col = find_header_index(header_row, ["-3sd"], start=hqc_start, end=hqc_end)
    hqc_cv_col = find_header_index(header_row, ["%cv"], start=hqc_start, end=hqc_end)

    lqc_mean_col = find_header_index(header_row, ["qc mean"], start=lqc_start, end=lqc_end)
    lqc_upper2_col = find_header_index(header_row, ["+2sd"], start=lqc_start, end=lqc_end)
    lqc_lower2_col = find_header_index(header_row, ["-2sd"], start=lqc_start, end=lqc_end)
    lqc_upper3_col = find_header_index(header_row, ["+3sd"], start=lqc_start, end=lqc_end)
    lqc_lower3_col = find_header_index(header_row, ["-3sd"], start=lqc_start, end=lqc_end)
    lqc_cv_col = find_header_index(header_row, ["%cv"], start=lqc_start, end=lqc_end)

    targets = []
    hqc_mean = value_row[hqc_mean_col] if hqc_mean_col is not None else None
    hqc_sd = calculate_sd(
        hqc_mean,
        upper2=value_row[hqc_upper2_col] if hqc_upper2_col is not None else None,
        lower2=value_row[hqc_lower2_col] if hqc_lower2_col is not None else None,
        upper3=value_row[hqc_upper3_col] if hqc_upper3_col is not None else None,
        lower3=value_row[hqc_lower3_col] if hqc_lower3_col is not None else None,
        cv=value_row[hqc_cv_col] if hqc_cv_col is not None else None,
    )
    if hqc_sd is not None and hqc_mean is not None and not pd.isna(hqc_mean):
        targets.append({
            "analyte": analyte,
            "qc_level": "High",
            "target_mean": float(hqc_mean),
            "target_sd": float(hqc_sd),
            "effective_from": file_date,
        })

    lqc_mean = value_row[lqc_mean_col] if lqc_mean_col is not None else None
    lqc_sd = calculate_sd(
        lqc_mean,
        upper2=value_row[lqc_upper2_col] if lqc_upper2_col is not None else None,
        lower2=value_row[lqc_lower2_col] if lqc_lower2_col is not None else None,
        upper3=value_row[lqc_upper3_col] if lqc_upper3_col is not None else None,
        lower3=value_row[lqc_lower3_col] if lqc_lower3_col is not None else None,
        cv=value_row[lqc_cv_col] if lqc_cv_col is not None else None,
    )
    if lqc_sd is not None and lqc_mean is not None and not pd.isna(lqc_mean):
        targets.append({
            "analyte": analyte,
            "qc_level": "Low",
            "target_mean": float(lqc_mean),
            "target_sd": float(lqc_sd),
            "effective_from": file_date,
        })

    return targets


def find_qc_run_header_row(df):
    for row_idx in range(min(len(df), 80)):
        row_vec = df.iloc[row_idx].tolist()
        row_strs = [str(v).strip().lower() if pd.notna(v) else "" for v in row_vec]
        if "run" in row_strs and "date" in row_strs and row_strs.count("result") >= 2:
            return row_idx, row_strs
    return None, None


def parse_qc_run_rows_from_sheet(sheet_name, df):
    analyte = find_analyte_name_in_workbook(sheet_name, sample_row=df.iloc[0].tolist() if len(df) > 0 else None)
    if analyte is None:
        return []

    header_row_idx, header_row = find_qc_run_header_row(df)
    if header_row_idx is None:
        return []

    date_col = find_header_index(header_row, ["date"])
    run_col = find_header_index(header_row, ["run"])
    result_cols = [i for i, value in enumerate(header_row) if value == "result"]
    if len(result_cols) < 2:
        return []

    hqc_result_col = result_cols[0]
    lqc_result_col = result_cols[1]

    records = []
    for row_idx in range(header_row_idx + 1, len(df)):
        row = df.iloc[row_idx]
        if pd.isna(row.iloc[hqc_result_col]) and pd.isna(row.iloc[lqc_result_col]):
            continue

        run_date = parse_date_value(row.iloc[date_col]) if date_col is not None else None
        if run_date is None:
            continue

        replicate = None
        if run_col is not None:
            try:
                replicate_val = row.iloc[run_col]
                if pd.notna(replicate_val):
                    replicate = int(float(replicate_val))
            except Exception:
                replicate = None

        for qc_level, col_idx in [("High", hqc_result_col), ("Low", lqc_result_col)]:
            raw_value = row.iloc[col_idx]
            if pd.isna(raw_value) or str(raw_value).strip() == "":
                continue
            try:
                concentration = float(raw_value)
            except Exception:
                continue

            records.append({
                "analyte": analyte,
                "qc_level": qc_level,
                "run_date": run_date,
                "concentration": concentration,
                "replicate": replicate if replicate is not None else 1,
            })

    return records


def get_qc_target(analyte_name, qc_level, as_of_date=None, db_path=None):
    db_path = db_path or DB_PATH
    if not Path(db_path).exists():
        return None

    # Normalize as_of_date to YYYY-MM-DD string
    if as_of_date is None:
        as_of_date = datetime.today().strftime("%Y-%m-%d")
    else:
        try:
            # Accept date/datetime objects or strings
            if not isinstance(as_of_date, str):
                as_of_date = pd.to_datetime(as_of_date, errors='coerce').strftime("%Y-%m-%d")
            else:
                parsed = pd.to_datetime(as_of_date, errors='coerce')
                if pd.isna(parsed):
                    as_of_date = datetime.today().strftime("%Y-%m-%d")
                else:
                    as_of_date = parsed.strftime("%Y-%m-%d")
        except Exception:
            as_of_date = datetime.today().strftime("%Y-%m-%d")

    conn = sqlite3.connect(str(db_path))
    # Primary lookup: effective range covering the as_of_date
    query = """
        SELECT qt.target_mean, qt.target_sd
        FROM qc_targets qt
        JOIN analytes a ON qt.analyte_id = a.analyte_id
        WHERE a.name = ?
          AND qt.qc_level = ?
          AND qt.effective_from <= ?
          AND (qt.effective_to IS NULL OR qt.effective_to >= ?)
        ORDER BY qt.effective_from DESC
        LIMIT 1
    """
    row = conn.execute(query, (analyte_name, qc_level, as_of_date, as_of_date)).fetchone()
    if not row:
        # Fallback 1: most recent effective_from <= as_of_date (ignore effective_to)
        query2 = """
            SELECT qt.target_mean, qt.target_sd
            FROM qc_targets qt
            JOIN analytes a ON qt.analyte_id = a.analyte_id
            WHERE a.name = ?
              AND qt.qc_level = ?
              AND qt.effective_from <= ?
            ORDER BY qt.effective_from DESC
            LIMIT 1
        """
        row = conn.execute(query2, (analyte_name, qc_level, as_of_date)).fetchone()

    if not row:
        # Fallback 2: most recent target regardless of date
        query3 = """
            SELECT qt.target_mean, qt.target_sd
            FROM qc_targets qt
            JOIN analytes a ON qt.analyte_id = a.analyte_id
            WHERE a.name = ?
              AND qt.qc_level = ?
            ORDER BY qt.effective_from DESC
            LIMIT 1
        """
        row = conn.execute(query3, (analyte_name, qc_level)).fetchone()

    conn.close()
    if row:
        return {"target_mean": row[0], "target_sd": row[1]}
    return None


def get_per_date_targets(analyte_name, qc_level, run_dates, db_path=None):
    """Return a list of (mean, sd) tuples — one per run_date — using the target active on each date."""
    return [
        get_qc_target(analyte_name, qc_level, as_of_date=d, db_path=db_path)
        for d in run_dates
    ]


def get_all_qc_targets(db_path=None):
    """Return all rows from qc_targets joined with analyte names."""
    db_path = db_path or DB_PATH
    if not Path(db_path).exists():
        return pd.DataFrame()
    conn = sqlite3.connect(str(db_path))
    df = pd.read_sql_query("""
        SELECT qt.target_id, a.name AS analyte, qt.qc_level, qt.lot_number,
               qt.target_mean, qt.target_sd, qt.effective_from, qt.effective_to
        FROM qc_targets qt
        JOIN analytes a ON qt.analyte_id = a.analyte_id
        ORDER BY a.name, qt.qc_level, qt.effective_from
    """, conn)
    conn.close()
    return df


def insert_qc_target(analyte_name, qc_level, target_mean, target_sd,
                     effective_from, effective_to=None, lot_number=None, db_path=None):
    """Insert or replace a QC target row. Also closes the previous open row for that analyte/level."""
    db_path = db_path or DB_PATH
    ensure_db_initialized(db_path)
    conn = get_connection(db_path)
    cursor = conn.cursor()

    cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (analyte_name,))
    row = cursor.fetchone()
    if not row:
        conn.close()
        raise ValueError(f"Analyte '{analyte_name}' not found in database.")
    analyte_id = row[0]

    # Close the previously open-ended row for this analyte+level (set effective_to = effective_from - 1 day)
    cursor.execute("""
        UPDATE qc_targets
        SET effective_to = date(?, '-1 day')
        WHERE analyte_id = ? AND qc_level = ? AND effective_to IS NULL
          AND effective_from < ?
    """, (effective_from, analyte_id, qc_level, effective_from))

    cursor.execute("""
        INSERT INTO qc_targets (analyte_id, qc_level, lot_number, target_mean, target_sd, effective_from, effective_to)
        VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(analyte_id, qc_level, lot_number, effective_from) DO UPDATE SET
            target_mean=excluded.target_mean,
            target_sd=excluded.target_sd,
            effective_to=excluded.effective_to
    """, (analyte_id, qc_level, lot_number or "", target_mean, target_sd, effective_from, effective_to))

    conn.commit()
    conn.close()


def import_qc_targets_file(file_bytes, filename, db_path=None):
    """Import mean/SD targets from CSV/Excel file into qc_targets."""
    db_path = db_path or DB_PATH
    ensure_db_initialized(db_path)

    suffix = Path(filename).suffix.lower()
    if suffix == ".csv":
        df = pd.read_csv(BytesIO(file_bytes))
    elif suffix in {".xls", ".xlsx"}:
        df = pd.read_excel(BytesIO(file_bytes))
    else:
        raise ValueError("Unsupported target file type. Please upload CSV or Excel.")

    if df.empty:
        return "No target rows found in file."

    col_map = {str(c).strip().lower(): c for c in df.columns}

    def pick(*names):
        for n in names:
            if n in col_map:
                return col_map[n]
        return None

    analyte_col = pick("analyte", "hormone", "compound", "name")
    level_col = pick("qc_level", "qc level", "level", "type")
    mean_col = pick("target_mean", "target mean", "mean", "qc mean")
    sd_col = pick("target_sd", "target sd", "sd")
    from_col = pick("effective_from", "effective from", "from", "date")
    to_col = pick("effective_to", "effective to", "to")
    lot_col = pick("lot_number", "lot number", "lot")

    if not analyte_col or not level_col or not mean_col or not sd_col:
        raise ValueError(
            "Targets file must include columns for analyte, qc_level, target_mean, and target_sd."
        )

    default_from = extract_date_from_filename(filename) or datetime.today().strftime("%Y-%m-%d")
    imported = 0
    skipped = 0

    for _, row in df.iterrows():
        analyte_name = str(row[analyte_col]).strip() if pd.notna(row[analyte_col]) else ""
        if not analyte_name:
            skipped += 1
            continue

        qc_level = normalize_qc_level(row[level_col])
        if qc_level is None:
            skipped += 1
            continue

        try:
            target_mean = float(row[mean_col])
            target_sd = float(row[sd_col])
        except Exception:
            skipped += 1
            continue

        effective_from = parse_date_value(row[from_col]) if from_col is not None else None
        effective_from = effective_from or default_from
        effective_to = parse_date_value(row[to_col]) if to_col is not None else None
        lot_number = str(row[lot_col]).strip() if lot_col is not None and pd.notna(row[lot_col]) else None

        insert_qc_target(
            analyte_name=analyte_name,
            qc_level=qc_level,
            target_mean=target_mean,
            target_sd=target_sd,
            effective_from=effective_from,
            effective_to=effective_to,
            lot_number=lot_number,
            db_path=db_path,
        )
        imported += 1

    return f"Imported/updated {imported} target rows" + (f" (skipped {skipped})" if skipped else "")


def import_excel_qc_file(file_bytes, filename, db_path=None, uploaded_by=None):
    """Import QC measurement or mean-value Excel data into the QC database."""
    ensure_db_initialized(db_path)
    conn = get_connection(db_path)
    cursor = conn.cursor()

    uploaded_by = str(uploaded_by).strip().upper() if uploaded_by else None

    source_filename = Path(filename).name
    cursor.execute("SELECT run_id FROM runs WHERE source_filename = ?", (source_filename,))
    if cursor.fetchone():
        conn.close()
        return f"Already imported: {source_filename}"

    excel_file = pd.ExcelFile(BytesIO(file_bytes))
    file_date = extract_date_from_filename(source_filename)

    analyte_sheets = []
    for sheet in excel_file.sheet_names:
        analyte = find_analyte_name_in_workbook(sheet)
        if analyte is not None:
            analyte_sheets.append((sheet, analyte))

    records = []
    qc_targets = []
    if analyte_sheets:
        sheet_date = file_date
        for sheet_name, analyte in analyte_sheets:
            df_sheet = pd.read_excel(excel_file, sheet_name=sheet_name, header=None)
            sheet_records = parse_qc_run_rows_from_sheet(sheet_name, df_sheet)
            records.extend(sheet_records)
            if sheet_date is None and sheet_records:
                found_dates = [r["run_date"] for r in sheet_records if r.get("run_date")]
                if found_dates:
                    sheet_date = min(found_dates)
        for sheet_name, analyte in analyte_sheets:
            df_sheet = pd.read_excel(excel_file, sheet_name=sheet_name, header=None)
            qc_targets.extend(parse_qc_targets_from_sheet(sheet_name, df_sheet, sheet_date or datetime.today().strftime("%Y-%m-%d")))
    else:
        df_raw = pd.read_excel(BytesIO(file_bytes))
        if df_raw.empty:
            conn.close()
            return "The Excel file contains no data."

        df = df_raw.copy()
        analyte_col = find_column(df, ["analyte", "compound", "name", "assay"])
        date_col = find_column(df, ["date", "run date", "run_date", "measurement date", "sample date"])
        qc_level_col = find_column(df, ["qc level", "level", "type", "sample type"])
        value_col = find_column(df, ["concentration", "mean", "value", "result", "measurement"])
        replicate_col = find_column(df, ["replicate", "rep", "replicate number"])

        hqc_value_columns = [col for col in df.columns if any(token in str(col).strip().lower() for token in ["hqc", "high"]) and any(token in str(col).strip().lower() for token in ["mean", "conc"])]
        lqc_value_columns = [col for col in df.columns if any(token in str(col).strip().lower() for token in ["lqc", "low"]) and any(token in str(col).strip().lower() for token in ["mean", "conc"])]

        if analyte_col is None:
            conn.close()
            raise ValueError("Excel file must include an analyte column (for example, Analyte, Compound, or Name).")

        if hqc_value_columns or lqc_value_columns:
            for _, row in df.iterrows():
                analyte = str(row[analyte_col]).strip() if pd.notna(row[analyte_col]) else None
                if not analyte:
                    continue

                run_date = parse_date_value(row[date_col]) if date_col is not None else None
                run_date = run_date or file_date
                if run_date is None:
                    raise ValueError("Excel file must include a date column or the filename must contain a date.")

                for col in hqc_value_columns:
                    concentration = row[col]
                    if pd.isna(concentration) or str(concentration).strip() == "":
                        continue
                    records.append({
                        "analyte": analyte,
                        "qc_level": "High",
                        "run_date": run_date,
                        "concentration": float(concentration),
                        "replicate": int(row[replicate_col]) if replicate_col and pd.notna(row[replicate_col]) else 1,
                    })

                for col in lqc_value_columns:
                    concentration = row[col]
                    if pd.isna(concentration) or str(concentration).strip() == "":
                        continue
                    records.append({
                        "analyte": analyte,
                        "qc_level": "Low",
                        "run_date": run_date,
                        "concentration": float(concentration),
                        "replicate": int(row[replicate_col]) if replicate_col and pd.notna(row[replicate_col]) else 1,
                    })
        elif qc_level_col is not None and value_col is not None:
            for _, row in df.iterrows():
                analyte = str(row[analyte_col]).strip() if pd.notna(row[analyte_col]) else None
                if not analyte:
                    continue

                qc_level = normalize_qc_level(row[qc_level_col])
                if qc_level is None:
                    continue

                run_date = parse_date_value(row[date_col]) if date_col is not None else None
                run_date = run_date or file_date
                if run_date is None:
                    raise ValueError("Excel file must include a date column or the filename must contain a date.")

                concentration = row[value_col]
                if pd.isna(concentration) or str(concentration).strip() == "":
                    continue

                records.append({
                    "analyte": analyte,
                    "qc_level": qc_level,
                    "run_date": run_date,
                    "concentration": float(concentration),
                    "replicate": int(row[replicate_col]) if replicate_col and pd.notna(row[replicate_col]) else 1,
                })
        else:
            conn.close()
            raise ValueError(
                "Excel import requires either HQC/LQC value columns (e.g. HQC Mean, LQC Mean) or a QC Level column plus a concentration/mean column."
            )

    if not records and not qc_targets:
        conn.close()
        return "No QC records were found in the Excel file."

    run_date = records[0]["run_date"] if records else file_date or datetime.today().strftime("%Y-%m-%d")
    cursor.execute(
        "INSERT INTO runs (run_date, panel, source_filename, method_name, data_path, uploaded_by) VALUES (?, ?, ?, ?, ?, ?)",
        (run_date, 1, source_filename, None, None, uploaded_by)
    )
    run_id = cursor.lastrowid

    if qc_targets:
        for target in qc_targets:
            cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (target["analyte"],))
            row = cursor.fetchone()
            if not row:
                continue
            analyte_id = row[0]
            cursor.execute(
                "INSERT OR REPLACE INTO qc_targets (analyte_id, qc_level, lot_number, target_mean, target_sd, effective_from) VALUES (?, ?, ?, ?, ?, ?)",
                (analyte_id, target["qc_level"], None, target["target_mean"], target["target_sd"], target["effective_from"])
            )

    analyte_id_map = {}
    for record in records:
        analyte = record["analyte"]
        cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (analyte,))
        row = cursor.fetchone()
        if row:
            analyte_id_map[analyte] = row[0]

    cursor.execute("SELECT type_code, type_id FROM sample_types")
    type_map = dict(cursor.fetchall())

    imported_count = 0
    for index, record in enumerate(records, start=1):
        analyte = record["analyte"]
        if analyte not in analyte_id_map:
            continue

        qc_level = record["qc_level"]
        concentration = record["concentration"]
        run_date = record["run_date"]
        replicate = record["replicate"]

        data_filename = f"QC_{qc_level}_{analyte}_{run_date}_{replicate}"
        data_filename = re.sub(r"[^A-Za-z0-9_.-]", "_", data_filename)

        sample_type_id = type_map["qc"]
        cursor.execute(
            "INSERT INTO samples (run_id, data_filename, sample_name, sample_type_id, collection_date, qc_level, qc_replicate) VALUES (?, ?, ?, ?, ?, ?, ?)",
            (run_id, data_filename, analyte, sample_type_id, run_date, qc_level, replicate)
        )
        sample_id = cursor.lastrowid

        cursor.execute(
            "INSERT OR REPLACE INTO results (sample_id, analyte_id, concentration) VALUES (?, ?, ?)",
            (sample_id, analyte_id_map[analyte], concentration)
        )
        imported_count += 1

    conn.commit()
    conn.close()
    return f"Imported {imported_count} QC records from {source_filename}."


# ==============================================================================
# QC DATA QUERIES
# ==============================================================================

def get_qc_data(db_path=None):
    """Pull all QC results from the database."""
    db_path = db_path or DB_PATH
    if not db_path.exists():
        return pd.DataFrame()

    conn = get_connection(db_path)
    query = """
        SELECT
            r.run_date,
            a.name as analyte,
            s.qc_level,
            r.uploaded_by,
            AVG(res.concentration) as concentration
        FROM results res
        JOIN samples s ON res.sample_id = s.sample_id
        JOIN runs r ON s.run_id = r.run_id
        JOIN analytes a ON res.analyte_id = a.analyte_id
        JOIN sample_types st ON s.sample_type_id = st.type_id
        WHERE st.type_code = 'qc'
          AND res.concentration IS NOT NULL
        GROUP BY r.run_date, a.name, s.qc_level, r.uploaded_by
        ORDER BY a.name, r.run_date
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


def query_run_summary(db_path=None):
    """Get summary of all runs."""
    db_path = db_path or DB_PATH
    if not db_path.exists():
        return pd.DataFrame()

    conn = get_connection(db_path)
    query = """
        SELECT
            r.run_id, r.run_date, r.panel, r.method_name, r.source_filename, r.uploaded_by, r.imported_at,
            COUNT(DISTINCT s.sample_id) as sample_count
        FROM runs r
        LEFT JOIN samples s ON r.run_id = s.run_id
        GROUP BY r.run_id ORDER BY r.run_date DESC
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


# ==============================================================================
# QC EXPORT FUNCTIONS
# ==============================================================================

def format_date(date_str):
    """Convert YYYY-MM-DD to DD/MM/YYYY."""
    parts = date_str.split("-")
    return f"{parts[2]}/{parts[1]}/{parts[0]}"


def export_hormone_csv(analyte_name, hqc_data, lqc_data):
    """Create CSV data for one hormone with HQC and LQC side by side."""
    hqc_mean = hqc_data["concentration"].mean() if not hqc_data.empty else np.nan
    hqc_sd = hqc_data["concentration"].std() if len(hqc_data) > 1 else np.nan

    lqc_mean = lqc_data["concentration"].mean() if not lqc_data.empty else np.nan
    lqc_sd = lqc_data["concentration"].std() if len(lqc_data) > 1 else np.nan

    all_dates = sorted(set(hqc_data["run_date"].tolist() + lqc_data["run_date"].tolist()))

    hqc_by_date = dict(zip(hqc_data["run_date"], hqc_data["concentration"]))
    lqc_by_date = dict(zip(lqc_data["run_date"], lqc_data["concentration"]))

    rows = []
    uploader_by_date = {}
    for row_info in pd.concat([hqc_data, lqc_data]).to_dict("records"):
        run_date = row_info.get("run_date")
        if run_date:
            uploader_by_date[run_date] = row_info.get("uploaded_by", "")

    for date in all_dates:
        row = {"Date": format_date(date), "User_Initials": uploader_by_date.get(date, "")}

        hqc_conc = hqc_by_date.get(date)
        row["HQC_Conc"] = round(hqc_conc, 4) if hqc_conc is not None else ""
        row["HQC_Mean"] = round(hqc_mean, 4) if not np.isnan(hqc_mean) else ""
        row["HQC_+2SD"] = round(hqc_mean + 2 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""
        row["HQC_-2SD"] = round(hqc_mean - 2 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""
        row["HQC_+3SD"] = round(hqc_mean + 3 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""
        row["HQC_-3SD"] = round(hqc_mean - 3 * hqc_sd, 4) if not np.isnan(hqc_sd) else ""

        lqc_conc = lqc_by_date.get(date)
        row["LQC_Conc"] = round(lqc_conc, 4) if lqc_conc is not None else ""
        row["LQC_Mean"] = round(lqc_mean, 4) if not np.isnan(lqc_mean) else ""
        row["LQC_+2SD"] = round(lqc_mean + 2 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""
        row["LQC_-2SD"] = round(lqc_mean - 2 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""
        row["LQC_+3SD"] = round(lqc_mean + 3 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""
        row["LQC_-3SD"] = round(lqc_mean - 3 * lqc_sd, 4) if not np.isnan(lqc_sd) else ""

        rows.append(row)

    return pd.DataFrame(rows)


# ==============================================================================
# QC CHART FUNCTIONS
# ==============================================================================

def flag_outliers(concentrations, sd2_upper, sd2_lower, sd3_upper, sd3_lower):
    """Flag values exceeding 2SD or 3SD bands."""
    flags = [False] * len(concentrations)

    for i, conc in enumerate(concentrations):
        if conc > sd3_upper or conc < sd3_lower:
            flags[i] = True

        if conc > sd2_upper or conc < sd2_lower:
            if i > 0 and (concentrations[i - 1] > sd2_upper or concentrations[i - 1] < sd2_lower):
                flags[i] = True
                flags[i - 1] = True

    return flags


def make_qc_chart(dates, concentrations, mean_val, sd2_upper, sd2_lower, sd3_upper, sd3_lower, title, uploader_initials=None, means_per_point=None, sds_per_point=None):
    """Create Levey-Jennings chart with 2SD/3SD bands.
    If means_per_point/sds_per_point are supplied they override the flat mean_val/sdX_upper/lower
    and the reference lines are drawn as step functions."""
    # Build per-point SD band values
    if means_per_point and sds_per_point and len(means_per_point) == len(dates):
        pp_mean    = [float(m) for m in means_per_point]
        pp_sd2_u   = [float(m) + 2 * float(s) for m, s in zip(means_per_point, sds_per_point)]
        pp_sd2_l   = [float(m) - 2 * float(s) for m, s in zip(means_per_point, sds_per_point)]
        pp_sd3_u   = [float(m) + 3 * float(s) for m, s in zip(means_per_point, sds_per_point)]
        pp_sd3_l   = [float(m) - 3 * float(s) for m, s in zip(means_per_point, sds_per_point)]
        # also recalculate the flat values used for flagging as the last active target
        mean_val   = pp_mean[-1]
        sd2_upper  = pp_sd2_u[-1]
        sd2_lower  = pp_sd2_l[-1]
        sd3_upper  = pp_sd3_u[-1]
        sd3_lower  = pp_sd3_l[-1]
    else:
        pp_mean  = [mean_val]  * len(dates)
        pp_sd2_u = [sd2_upper] * len(dates)
        pp_sd2_l = [sd2_lower] * len(dates)
        pp_sd3_u = [sd3_upper] * len(dates)
        pp_sd3_l = [sd3_lower] * len(dates)

    flags = flag_outliers(concentrations, sd2_upper, sd2_lower, sd3_upper, sd3_lower)
    initials_list = [(str(x).strip().upper() if pd.notna(x) and str(x).strip() else "NA") for x in (uploader_initials if uploader_initials is not None else [None] * len(dates))]
    if len(initials_list) != len(dates):
        initials_list = [initials_list[0] if initials_list else "NA"] * len(dates)

    fig = go.Figure()

    # Step reference lines (change at each effective date)
    fig.add_trace(go.Scatter(
        x=dates, y=pp_mean,
        mode="lines",
        line=dict(color="#008000", dash="dash", width=2, shape="hv"),
        name="Mean",
        hoverinfo="skip",
    ))
    fig.add_trace(go.Scatter(
        x=dates, y=pp_sd2_u,
        mode="lines",
        line=dict(color="#ff9800", dash="dot", width=1, shape="hv"),
        name="+2SD",
        hoverinfo="skip",
    ))
    fig.add_trace(go.Scatter(
        x=dates, y=pp_sd2_l,
        mode="lines",
        line=dict(color="#ff9800", dash="dot", width=1, shape="hv"),
        name="-2SD",
        hoverinfo="skip",
    ))
    fig.add_trace(go.Scatter(
        x=dates, y=pp_sd3_u,
        mode="lines",
        line=dict(color="#ff3d00", dash="dash", width=1, shape="hv"),
        name="+3SD",
        hoverinfo="skip",
    ))
    fig.add_trace(go.Scatter(
        x=dates, y=pp_sd3_l,
        mode="lines",
        line=dict(color="#ff3d00", dash="dash", width=1, shape="hv"),
        name="-3SD",
        hoverinfo="skip",
    ))

    # Per-point hover data includes that point's active mean/SD
    customdata_main = [
        [ini, pm, psu, psl, pu3, pl3]
        for ini, pm, psu, psl, pu3, pl3
        in zip(initials_list, pp_mean, pp_sd2_u, pp_sd2_l, pp_sd3_u, pp_sd3_l)
    ]
    fig.add_trace(go.Scatter(
        x=dates,
        y=concentrations,
        mode="lines+markers",
        marker=dict(size=9, color="#1976d2"),
        line=dict(color="#1976d2", width=2),
        name="Concentration",
        customdata=customdata_main,
        hovertemplate=(
            "Date: %{x}<br>"
            "Concentration: %{y:.3f}<br>"
            "Initials: %{customdata[0]}<br>"
            "Mean: %{customdata[1]:.3f}<br>"
            "+2SD: %{customdata[2]:.3f}<br>"
            "-2SD: %{customdata[3]:.3f}<br>"
            "+3SD: %{customdata[4]:.3f}<br>"
            "-3SD: %{customdata[5]:.3f}<extra></extra>"
        ),
    ))

    flagged_dates = [d for d, f in zip(dates, flags) if f]
    flagged_concs = [c for c, f in zip(concentrations, flags) if f]
    flagged_initials = [ini for ini, f in zip(initials_list, flags) if f]
    flagged_custom = [cd for cd, f in zip(customdata_main, flags) if f]

    if flagged_dates:
        fig.add_trace(go.Scatter(
            x=flagged_dates, y=flagged_concs,
            mode="markers",
            marker=dict(size=12, color="#d32f2f", symbol="triangle-up", line=dict(width=1, color="#b71c1c")),
            name="Flagged",
            customdata=flagged_custom,
            hovertemplate=(
                "Date: %{x}<br>"
                "Concentration: %{y:.3f}<br>"
                "Initials: %{customdata[0]}<br>"
                "Mean: %{customdata[1]:.3f}<br>"
                "+2SD: %{customdata[2]:.3f}<br>"
                "-2SD: %{customdata[3]:.3f}<br>"
                "+3SD: %{customdata[4]:.3f}<br>"
                "-3SD: %{customdata[5]:.3f}<extra></extra>"
            ),
        ))

    fig.update_layout(
        title=dict(text=title, y=0.97, yanchor="top"),
        xaxis_title="Date",
        yaxis_title="Concentration",
        height=420,
        plot_bgcolor="rgba(0,0,0,0)",
        paper_bgcolor="rgba(0,0,0,0)",
        font=dict(color="#666666"),
        legend=dict(orientation="h", yanchor="bottom", y=1.08, xanchor="right", x=1),
        margin=dict(t=115, b=40, l=60, r=20),
    )

    fig.update_xaxes(
        showgrid=True,
        gridcolor="rgba(200,200,200,0.3)",
        zeroline=False,
        tickangle=-45,
        title_standoff=10,
    )
    fig.update_yaxes(
        showgrid=True,
        gridcolor="rgba(200,200,200,0.3)",
        zeroline=False,
        tickformat=".3f",
        title_standoff=10,
    )

    return fig


def create_value_pictogram(concentrations, mean_val, sd_val):
    """Create an inline SVG pictogram (SD lines, history dots, recent circle)."""
    if not concentrations:
        svg = (
            "<svg xmlns='http://www.w3.org/2000/svg' width='220' height='56'>"
            "<rect x='1' y='1' width='218' height='54' rx='5' fill='none' stroke='#bdbdbd'/>"
            "<text x='110' y='34' text-anchor='middle' font-size='12' fill='#9e9e9e'>N/A</text>"
            "</svg>"
        )
        return f"data:image/svg+xml;utf8,{quote(svg)}"

    w = 220
    h = 56
    left = 10
    right = w - 10
    y_mid = 28
    y_hist = 40

    if sd_val and sd_val > 0 and not pd.isna(sd_val):
        effective_sd = float(sd_val)
    else:
        data_min = min(concentrations)
        data_max = max(concentrations)
        span = max(data_max - data_min, 1e-6)
        effective_sd = max(span / 6.0, abs(float(mean_val)) * 0.05, 1e-6)

    lo = float(mean_val) - 3 * effective_sd
    hi = float(mean_val) + 3 * effective_sd

    def x_of(value):
        ratio = (float(value) - lo) / max(hi - lo, 1e-9)
        ratio = min(max(ratio, 0.0), 1.0)
        return left + ratio * (right - left)

    x_mean = x_of(mean_val)
    x_recent = x_of(concentrations[-1])
    x_2sd_l = x_of(mean_val - 2 * effective_sd)
    x_2sd_r = x_of(mean_val + 2 * effective_sd)
    x_3sd_l = x_of(mean_val - 3 * effective_sd)
    x_3sd_r = x_of(mean_val + 3 * effective_sd)

    hist_vals = concentrations[:-1][-8:]
    hist_circles = "".join(
        f"<circle cx='{x_of(v):.1f}' cy='{y_hist}' r='3.8' fill='#455a64' opacity='0.9'/>"
        for v in hist_vals
    )

    svg = f"""
<svg xmlns='http://www.w3.org/2000/svg' width='{w}' height='{h}'>
  <rect x='1' y='1' width='{w - 2}' height='{h - 2}' rx='6' fill='none' stroke='#90a4ae' stroke-width='1'/>
  <line x1='{left}' y1='{y_mid}' x2='{right}' y2='{y_mid}' stroke='#90a4ae' stroke-width='1.2'/>
  <line x1='{x_3sd_l:.1f}' y1='14' x2='{x_3sd_l:.1f}' y2='42' stroke='#ef5350' stroke-width='2'/>
  <line x1='{x_3sd_r:.1f}' y1='14' x2='{x_3sd_r:.1f}' y2='42' stroke='#ef5350' stroke-width='2'/>
  <line x1='{x_2sd_l:.1f}' y1='17' x2='{x_2sd_l:.1f}' y2='39' stroke='#ff9800' stroke-width='1.8'/>
  <line x1='{x_2sd_r:.1f}' y1='17' x2='{x_2sd_r:.1f}' y2='39' stroke='#ff9800' stroke-width='1.8'/>
  <line x1='{x_mean:.1f}' y1='12' x2='{x_mean:.1f}' y2='44' stroke='#111111' stroke-width='2.2'/>
  {hist_circles}
  <circle cx='{x_recent:.1f}' cy='{y_mid}' r='8.5' fill='#fff200' stroke='#111111' stroke-width='1.6'/>
  <circle cx='{x_mean:.1f}' cy='{y_mid}' r='3.3' fill='#111111'/>
</svg>
""".strip()

    return f"data:image/svg+xml;utf8,{quote(svg)}"


def generate_final_report(db_path=None):
    """Generate comprehensive QC report for all hormones."""
    db_path = db_path or DB_PATH
    if not db_path.exists():
        return None
    
    df_qc = get_qc_data(db_path)
    if df_qc.empty:
        return None
    
    analytes = sorted(df_qc["analyte"].unique())
    report_data = []
    
    for analyte in analytes:
        analyte_data = df_qc[df_qc["analyte"] == analyte]
        dashboard_url = f"?mode=Dashboard&hormone={quote(str(analyte))}"
        
        for qc_level in ["High", "Low"]:
            level_data = analyte_data[analyte_data["qc_level"] == qc_level].reset_index(drop=True)
            
            if level_data.empty:
                report_data.append({
                    "Hormone": analyte,
                    "Go to Dashboard": dashboard_url,
                    "QC Level": "HQC" if qc_level == "High" else "LQC",
                    "Pictogram": create_value_pictogram([], 0.0, 0.0),
                    "Recent": "—",
                    "Mean": "—",
                    "Min": "—",
                    "Max": "—",
                    "N": 0,
                    "Status": "N/A"
                })
                continue
            
            concentrations = level_data["concentration"].tolist()
            recent = concentrations[-1]
            mean_val = level_data["concentration"].mean()
            sd = level_data["concentration"].std() if len(level_data) > 1 else 0
            
            target = get_qc_target(analyte, qc_level, as_of_date=level_data["run_date"].max())
            if target:
                mean_val = float(target["target_mean"])
                sd = float(target["target_sd"])
            
            sd2_upper = mean_val + 2 * sd
            sd2_lower = mean_val - 2 * sd

            if pd.isna(sd) or float(sd) <= 0:
                status = "N/A (SD unavailable)"
            else:
                z_abs = abs((recent - mean_val) / sd)
                if z_abs > 3:
                    status = "⚠ Out of Range (>3 SD)"
                elif z_abs > 2:
                    status = "⚠ Out of Range (>2 SD)"
                else:
                    status = "✓ OK (within 2 SD)"
            
            report_data.append({
                "Hormone": analyte,
                "Go to Dashboard": dashboard_url,
                "QC Level": "HQC" if qc_level == "High" else "LQC",
                "Pictogram": create_value_pictogram(concentrations, mean_val, sd),
                "Recent": f"{recent:.3f}",
                "Mean": f"{mean_val:.3f}",
                "Min": f"{min(concentrations):.3f}",
                "Max": f"{max(concentrations):.3f}",
                "N": len(concentrations),
                "Status": status
            })
    
    return pd.DataFrame(report_data)


# ==============================================================================
# STREAMLIT APP
# ==============================================================================

def main():
    st.set_page_config(page_title="QC Studio", layout="wide")
    st.title("🧪 QC Studio")
    st.markdown("Integrated steroid panel database, QC export, and dashboard platform")

    # Sidebar navigation
    module_options = ["Dashboard", "Database", "Export", "Report"]
    query_mode = str(st.query_params.get("mode", "")).strip()
    default_mode_idx = module_options.index(query_mode) if query_mode in module_options else 0
    mode = st.sidebar.radio(
        "Select Module",
        module_options,
        index=default_mode_idx,
        help="Choose between viewing QC charts, managing the database, exporting data, or generating a final report"
    )
    st.query_params["mode"] = mode

    st.sidebar.markdown("---")
    st.sidebar.markdown(f"**Database:** `{DB_PATH.name}`")

    if mode == "Database":
        st.header("📊 Steroid Panel Database")

        st.info("💡 **How it works:** Upload CSV or Excel files to automatically create and populate the database. The schema and reference data are generated on-demand from your first file upload.")

        col1, col2 = st.columns(2)

        with col1:
            initials = st.text_input(
                "Enter your initials",
                max_chars=6,
                key="qc_uploader_initials",
                help="Enter the initials of the user uploading this QC file. Leave blank if not available."
            )
            uploaded_file = st.file_uploader(
                "📁 Import QC Data (CSV or Excel)",
                type=["csv", "xls", "xlsx"],
                key="qc_data_uploader"
            )
            if uploaded_file:
                tmp_path = None
                try:
                    initials = str(initials).strip().upper() if initials else None
                    with st.spinner("Importing QC data..."):
                        suffix = Path(uploaded_file.name).suffix.lower()
                        if suffix == ".csv":
                            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as tmp:
                                tmp.write(uploaded_file.getbuffer())
                                tmp_path = tmp.name
                            result = import_csv(tmp_path, uploaded_by=initials)
                        elif suffix in {".xls", ".xlsx"}:
                            result = import_excel_qc_file(uploaded_file.read(), uploaded_file.name, uploaded_by=initials)
                        else:
                            raise ValueError("Unsupported file type. Upload a CSV or Excel file.")
                    st.success(result)
                    get_qc_data.clear() if hasattr(get_qc_data, 'clear') else None
                except Exception as e:
                    st.error(f"Import failed: {e}")
                finally:
                    if tmp_path and Path(tmp_path).exists():
                        os.unlink(tmp_path)

        st.markdown("---")

        st.subheader("📋 Run Summary")
        df_runs = query_run_summary()
        if df_runs.empty:
            st.info("No data imported yet. Upload a CSV file to get started.")
        else:
            st.dataframe(df_runs, use_container_width=True)

        st.markdown("---")
        st.subheader("🎯 QC Targets Manager")
        st.markdown("View existing mean/SD targets per hormone and add new ones when lot changes.")

        df_targets = get_all_qc_targets()
        if df_targets.empty:
            st.info("No QC targets stored yet. Import an Excel workbook or add one manually below.")
        else:
            st.dataframe(df_targets, use_container_width=True, hide_index=True)

        with st.expander("📤 Upload Mean/SD Targets File"):
            st.caption(
                "Upload CSV/Excel with columns: analyte, qc_level, target_mean, target_sd, "
                "and optional effective_from, effective_to, lot_number."
            )
            targets_file = st.file_uploader(
                "Upload targets file",
                type=["csv", "xls", "xlsx"],
                key="qc_targets_file_uploader",
            )
            if targets_file and st.button("⬆️ Import Targets File", use_container_width=True, key="import_targets_file_btn"):
                try:
                    msg = import_qc_targets_file(targets_file.read(), targets_file.name)
                    st.success(msg)
                    st.rerun()
                except Exception as e:
                    st.error(f"Targets import failed: {e}")

        with st.expander("➕ Add / Update QC Target"):
            if not DB_PATH.exists():
                st.warning("Import data first to initialise the database.")
            else:
                all_analyte_names = sorted([a["name"] for a in ANALYTES])
                t_col1, t_col2 = st.columns(2)
                with t_col1:
                    t_analyte = st.selectbox("Hormone", all_analyte_names, key="t_analyte")
                    t_level   = st.selectbox("QC Level", ["High (HQC)", "Low (LQC)"], key="t_level")
                    t_lot     = st.text_input("Lot Number (optional)", key="t_lot")
                with t_col2:
                    t_mean    = st.number_input("Target Mean", min_value=0.0, format="%.4f", key="t_mean")
                    t_sd      = st.number_input("Target SD",   min_value=0.0, format="%.4f", key="t_sd")
                    t_from    = st.date_input("Effective From", key="t_from")
                    t_to      = st.date_input("Effective To (leave blank = open-ended)", value=None, key="t_to")

                if st.button("💾 Save QC Target", use_container_width=True):
                    try:
                        level_code = "High" if "High" in t_level else "Low"
                        insert_qc_target(
                            analyte_name  = t_analyte,
                            qc_level      = level_code,
                            target_mean   = t_mean,
                            target_sd     = t_sd,
                            effective_from= str(t_from),
                            effective_to  = str(t_to) if t_to else None,
                            lot_number    = t_lot or None,
                        )
                        st.success(f"Saved target for {t_analyte} {level_code} effective {t_from}.")
                        st.rerun()
                    except Exception as e:
                        st.error(f"Failed to save: {e}")

    elif mode == "Export":
        st.header("📤 QC Export")

        if not DB_PATH.exists():
            st.error("Database not found. Import data first in the Database tab.")
            return

        df_qc = get_qc_data()
        if df_qc.empty:
            st.info("No QC data found in the database.")
            return

        analytes = sorted(df_qc["analyte"].unique())

        st.subheader("Export Hormone CSVs")
        st.markdown("Generate CSV files with HQC and LQC values for all hormones.")

        if st.button("📥 Generate All CSV Files", use_container_width=True):
            exported = []
            temp_dir = tempfile.mkdtemp()

            with st.spinner("Generating export files..."):
                for analyte in analytes:
                    analyte_data = df_qc[df_qc["analyte"] == analyte]
                    hqc_data = analyte_data[analyte_data["qc_level"] == "High"].reset_index(drop=True)
                    lqc_data = analyte_data[analyte_data["qc_level"] == "Low"].reset_index(drop=True)

                    if hqc_data.empty and lqc_data.empty:
                        continue

                    df_export = export_hormone_csv(analyte, hqc_data, lqc_data)
                    export_path = Path(temp_dir) / f"{analyte}_QC.csv"
                    df_export.to_csv(export_path, index=False)
                    exported.append((analyte, export_path))

            st.success(f"Generated {len(exported)} CSV files!")

            for analyte, export_path in exported:
                with open(export_path, "rb") as file:
                    st.download_button(
                        label=f"📥 Download {analyte}_QC.csv",
                        data=file.read(),
                        file_name=f"{analyte}_QC.csv",
                        mime="text/csv",
                        use_container_width=True
                    )

    elif mode == "Dashboard":
        st.header("📈 QC Dashboard")

        if not DB_PATH.exists():
            st.error("Database not found. Import data first in the Database tab.")
            return

        df = get_qc_data()
        if df.empty:
            st.warning("No QC data found in the database. Import data to view charts.")
            return

        analytes = sorted(df["analyte"].unique())
        query_hormone = str(st.query_params.get("hormone", "")).strip()
        default_hormone_idx = analytes.index(query_hormone) if query_hormone in analytes else 0
        selected = st.sidebar.radio("Select Hormone", analytes, index=default_hormone_idx)
        st.query_params["hormone"] = selected

        analyte_data = df[df["analyte"] == selected]
        hqc_data = analyte_data[analyte_data["qc_level"] == "High"].reset_index(drop=True)
        lqc_data = analyte_data[analyte_data["qc_level"] == "Low"].reset_index(drop=True)

        chart_cols = st.columns(2)

        if hqc_data.empty:
            chart_cols[0].info(f"No HQC data for {selected}.")
        else:
            hqc_concentrations = hqc_data["concentration"].tolist()
            hqc_raw_dates = hqc_data["run_date"].tolist()
            hqc_dates = [d.replace("-", "/") for d in hqc_raw_dates]
            hqc_targets = get_per_date_targets(selected, "High", hqc_raw_dates)
            has_targets = any(t is not None for t in hqc_targets)
            if has_targets:
                hqc_means = [float(t["target_mean"]) if t else hqc_data["concentration"].mean() for t in hqc_targets]
                hqc_sds   = [float(t["target_sd"])   if t else hqc_data["concentration"].std()  for t in hqc_targets]
                hqc_mean_val = hqc_means[-1]; hqc_sd = hqc_sds[-1]
                chart_cols[0].caption("Mean/SD step-lines follow active QC targets per run date.")
            else:
                hqc_mean_val = hqc_data["concentration"].mean()
                hqc_sd = hqc_data["concentration"].std()
                hqc_means = None; hqc_sds = None

            if pd.isna(hqc_sd) or hqc_sd == 0:
                chart_cols[0].warning(f"Not enough HQC data points for {selected} to compute SD.")
            else:
                hqc_fig = make_qc_chart(
                    hqc_dates, hqc_concentrations,
                    hqc_mean_val,
                    hqc_mean_val + 2 * hqc_sd,
                    hqc_mean_val - 2 * hqc_sd,
                    hqc_mean_val + 3 * hqc_sd,
                    hqc_mean_val - 3 * hqc_sd,
                    title=f"{selected} — HQC",
                    uploader_initials=hqc_data["uploaded_by"].fillna("NA").tolist(),
                    means_per_point=hqc_means,
                    sds_per_point=hqc_sds,
                )
                chart_cols[0].plotly_chart(hqc_fig, use_container_width=True)

        if lqc_data.empty:
            chart_cols[1].info(f"No LQC data for {selected}.")
        else:
            lqc_concentrations = lqc_data["concentration"].tolist()
            lqc_raw_dates = lqc_data["run_date"].tolist()
            lqc_dates = [d.replace("-", "/") for d in lqc_raw_dates]
            lqc_targets = get_per_date_targets(selected, "Low", lqc_raw_dates)
            has_lqc_targets = any(t is not None for t in lqc_targets)
            if has_lqc_targets:
                lqc_means = [float(t["target_mean"]) if t else lqc_data["concentration"].mean() for t in lqc_targets]
                lqc_sds   = [float(t["target_sd"])   if t else lqc_data["concentration"].std()  for t in lqc_targets]
                lqc_mean_val = lqc_means[-1]; lqc_sd = lqc_sds[-1]
                chart_cols[1].caption("Mean/SD step-lines follow active QC targets per run date.")
            else:
                lqc_mean_val = lqc_data["concentration"].mean()
                lqc_sd = lqc_data["concentration"].std()
                lqc_means = None; lqc_sds = None

            if pd.isna(lqc_sd) or lqc_sd == 0:
                chart_cols[1].warning(f"Not enough LQC data points for {selected} to compute SD.")
            else:
                lqc_fig = make_qc_chart(
                    lqc_dates, lqc_concentrations,
                    lqc_mean_val,
                    lqc_mean_val + 2 * lqc_sd,
                    lqc_mean_val - 2 * lqc_sd,
                    lqc_mean_val + 3 * lqc_sd,
                    lqc_mean_val - 3 * lqc_sd,
                    title=f"{selected} — LQC",
                    uploader_initials=lqc_data["uploaded_by"].fillna("NA").tolist(),
                    means_per_point=lqc_means,
                    sds_per_point=lqc_sds,
                )
                chart_cols[1].plotly_chart(lqc_fig, use_container_width=True)

        st.markdown("---")
        st.markdown("**HQC and LQC Statistics**")
        stats_cols = st.columns(2)

        if not hqc_data.empty:
            with stats_cols[0]:
                st.subheader("HQC Statistics")
                st.metric("HQC Mean", f"{hqc_data['concentration'].mean():.4f}")
                st.metric("HQC SD", f"{hqc_data['concentration'].std():.4f}")
                st.metric("HQC Min", f"{hqc_data['concentration'].min():.4f}")
                st.metric("HQC Max", f"{hqc_data['concentration'].max():.4f}")

        if not lqc_data.empty:
            with stats_cols[1]:
                st.subheader("LQC Statistics")
                st.metric("LQC Mean", f"{lqc_data['concentration'].mean():.4f}")
                st.metric("LQC SD", f"{lqc_data['concentration'].std():.4f}")
                st.metric("LQC Min", f"{lqc_data['concentration'].min():.4f}")
                st.metric("LQC Max", f"{lqc_data['concentration'].max():.4f}")

        st.caption("Select a different hormone from the sidebar list.")

    elif mode == "Report":
        st.header("📋 QC Final Report")
        
        if not DB_PATH.exists():
            st.error("Database not found. Import data first in the Database tab.")
            return
        
        report_df = generate_final_report()
        if report_df is None or report_df.empty:
            st.info("No QC data found in the database.")
            return
        
        st.markdown("**Summary of all hormones with LQC and HQC values**")
        
        # Display as a table with styling
        st.dataframe(
            report_df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Hormone": st.column_config.TextColumn("Hormone", width="medium"),
                "Go to Dashboard": st.column_config.LinkColumn("Open Dashboard", width="small", display_text="Open"),
                "QC Level": st.column_config.TextColumn("QC Level", width="small"),
                "Pictogram": st.column_config.ImageColumn("Pictogram", width="medium"),
                "Recent": st.column_config.TextColumn("Recent Value", width="small"),
                "Mean": st.column_config.TextColumn("Mean", width="small"),
                "Min": st.column_config.TextColumn("Min", width="small"),
                "Max": st.column_config.TextColumn("Max", width="small"),
                "N": st.column_config.NumberColumn("N", width="tiny"),
                "Status": st.column_config.TextColumn("Status", width="medium"),
            }
        )
        
        st.markdown("---")
        st.markdown("**Export Report**")
        
        if st.button("📥 Export Report as CSV", use_container_width=True):
            export_df = report_df.drop(columns=["Pictogram"], errors="ignore")
            csv = export_df.to_csv(index=False)
            st.download_button(
                label="Download Report CSV",
                data=csv,
                file_name=f"QC_Report_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv",
                mime="text/csv",
                use_container_width=True
            )


if __name__ == "__main__":
    main()
