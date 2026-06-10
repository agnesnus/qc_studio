"""
Steroid Panel LIMS - Database Setup and CSV Importer
=====================================================
SQLite database for LC-MS/MS steroid panel results.
Supports two CSV export formats:
  - Old: 20260206_Panel1_conc(in).csv (simple, concentrations only)
  - New: ISOSP-23_20260420_P1(Sheet1).csv (rich metadata + concentrations)

Usage:
    python steroid_panel_lims.py init
    python steroid_panel_lims.py import <csv_path>
    python steroid_panel_lims.py runs
    python steroid_panel_lims.py qc <analyte_name> <Low|High>
    python steroid_panel_lims.py eqa <scheme> <year> <round>
"""

import sqlite3
import re
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Optional

import pandas as pd

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

    # Calibrators: 'Cal 0', 'Cal A' ... 'Cal F'
    cal_match = re.match(r"^Cal\s+([0A-F])$", base)
    if cal_match:
        return SampleInfo(
            data_filename=data_filename,
            sample_type="calibrator",
            calibrator_level=cal_match.group(1),
        )

    # QC samples: 'QC_Low1', 'QC_High2'
    qc_match = re.match(r"^QC_(Low|High)(\d+)$", base)
    if qc_match:
        return SampleInfo(
            data_filename=data_filename,
            sample_type="qc",
            qc_level=qc_match.group(1),
            qc_replicate=int(qc_match.group(2)),
        )

    # Blanks: 'Blank1', 'Blank2', 'Blank3'
    if re.match(r"^Blank\d*$", base):
        return SampleInfo(data_filename=data_filename, sample_type="blank")

    # Process Blanks: 'PBlank1', 'PBlank2', 'PB1'
    if re.match(r"^(PBlank|PB)\d*$", base):
        return SampleInfo(data_filename=data_filename, sample_type="process_blank")

    # EQA old format: 'SKML2026_1A-a' (dash-separated replicate)
    eqa_match = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])-([a-z])$", base)
    if eqa_match:
        return SampleInfo(
            data_filename=data_filename,
            sample_type="eqa",
            eqa_scheme=eqa_match.group(1),
            eqa_year=int(eqa_match.group(2)),
            eqa_round=int(eqa_match.group(3)),
            eqa_sample_code=eqa_match.group(4),
            eqa_replicate=eqa_match.group(5),
        )

    # EQA new format with special suffix: 'SKML2026_4B_nowash2'
    eqa_special = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])_(.+)$", base)
    if eqa_special:
        return SampleInfo(
            data_filename=data_filename,
            sample_type="eqa",
            eqa_scheme=eqa_special.group(1),
            eqa_year=int(eqa_special.group(2)),
            eqa_round=int(eqa_special.group(3)),
            eqa_sample_code=eqa_special.group(4),
            eqa_replicate=eqa_special.group(5),
        )

    # EQA new format with replicate: 'SKML2026_3Ab' (lowercase appended)
    eqa_new_rep = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])([a-z])$", base)
    if eqa_new_rep:
        return SampleInfo(
            data_filename=data_filename,
            sample_type="eqa",
            eqa_scheme=eqa_new_rep.group(1),
            eqa_year=int(eqa_new_rep.group(2)),
            eqa_round=int(eqa_new_rep.group(3)),
            eqa_sample_code=eqa_new_rep.group(4),
            eqa_replicate=eqa_new_rep.group(5),
        )

    # EQA new format first injection: 'SKML2026_3A' (no replicate)
    eqa_new = re.match(r"^([A-Za-z]+)(\d{4})_(\d)([A-Z])$", base)
    if eqa_new:
        return SampleInfo(
            data_filename=data_filename,
            sample_type="eqa",
            eqa_scheme=eqa_new.group(1),
            eqa_year=int(eqa_new.group(2)),
            eqa_round=int(eqa_new.group(3)),
            eqa_sample_code=eqa_new.group(4),
            eqa_replicate=None,
        )

    # Patient samples: '20260122_SST01'
    patient_match = re.match(r"^(\d{8})_(\w+)$", base)
    if patient_match:
        date_str = patient_match.group(1)
        formatted_date = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
        return SampleInfo(
            data_filename=data_filename,
            sample_type="patient",
            collection_date=formatted_date,
            patient_sequence=patient_match.group(2),
        )

    return SampleInfo(data_filename=data_filename, sample_type="patient")


def classify_from_instrument_type(instrument_type: str, level: str, data_filename: str) -> SampleInfo:
    """Use instrument-assigned Type column to classify, falling back to regex for Sample type."""
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
            data_filename=data_filename,
            sample_type="qc",
            qc_level=qc_level,
            qc_replicate=qc_replicate,
        )

    # "Sample" type — use regex classifier to distinguish EQA from patient
    return classify_sample(data_filename)


# ==============================================================================
# DATABASE OPERATIONS
# ==============================================================================


def get_connection(db_path=None):
    db_path = db_path or DB_PATH
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys = ON")
    return conn


def init_db(db_path=None):
    """Create schema and seed reference data."""
    conn = get_connection(db_path)
    conn.executescript(SCHEMA_SQL)

    for analyte in ANALYTES:
        conn.execute(
            "INSERT OR IGNORE INTO analytes (name, panel, display_order) VALUES (?, ?, ?)",
            (analyte["name"], analyte["panel"], analyte["display_order"]),
        )

    for st in SAMPLE_TYPES:
        conn.execute(
            "INSERT OR IGNORE INTO sample_types (type_code, description) VALUES (?, ?)",
            (st["type_code"], st["description"]),
        )

    conn.commit()
    conn.close()
    print(f"Database initialized: {db_path or DB_PATH}")


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
    """Parse 'ISOSP-23_20260420_P1(Sheet1).csv' or 'ISOSP-23_20260420_P2-Neg(Sheet1).csv' → run metadata."""
    basename = Path(filepath).name
    match = re.match(r"^(.+?)_(\d{8})_P(\d+)(?:-\w+)?\(Sheet\d+\)\.csv$", basename)
    if not match:
        raise ValueError(f"Filename does not match new pattern: {basename}")
    method_name = match.group(1)
    date_str = match.group(2)
    run_date = f"{date_str[:4]}-{date_str[4:6]}-{date_str[6:8]}"
    panel = int(match.group(3))
    return {"run_date": run_date, "panel": panel, "source_filename": basename, "method_name": method_name}


def parse_filename(filepath: str, fmt: str) -> dict:
    """Route to appropriate filename parser."""
    if fmt == "old":
        return parse_filename_old(filepath)
    return parse_filename_new(filepath)


# ==============================================================================
# CSV IMPORTER — OLD FORMAT
# ==============================================================================


def import_csv_old(csv_path: str, db_path=None):
    """Import old-format CSV (simple header with Sample + analyte Results columns)."""
    conn = get_connection(db_path)
    cursor = conn.cursor()

    meta = parse_filename_old(csv_path)

    cursor.execute("SELECT run_id FROM runs WHERE source_filename = ?", (meta["source_filename"],))
    if cursor.fetchone():
        print(f"Already imported: {meta['source_filename']}")
        conn.close()
        return

    cursor.execute(
        "INSERT INTO runs (run_date, panel, source_filename, method_name) VALUES (?, ?, ?, ?)",
        (meta["run_date"], meta["panel"], meta["source_filename"], meta["method_name"]),
    )
    run_id = cursor.lastrowid

    df = pd.read_csv(csv_path, header=0)
    df = df.iloc[1:]  # Drop "Data File / Calc. Conc." row
    df = df.reset_index(drop=True)

    analyte_columns = []
    for col in df.columns[1:]:
        analyte_name = col.replace(" Results", "").strip()
        analyte_columns.append(analyte_name)

    analyte_id_map = {}
    for name in analyte_columns:
        cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (name,))
        row = cursor.fetchone()
        if row:
            analyte_id_map[name] = row[0]
        else:
            raise ValueError(f"Unknown analyte: {name}. Run 'init' first.")

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
            (
                run_id, data_filename, None, sample_type_id, None,
                None, None, None,
                info.collection_date, info.patient_sequence,
                info.calibrator_level, info.qc_level, info.qc_replicate,
                info.eqa_scheme, info.eqa_year, info.eqa_round,
                info.eqa_sample_code, info.eqa_replicate,
            ),
        )
        sample_id = cursor.lastrowid

        for i, analyte_name in enumerate(analyte_columns):
            raw_value = row.iloc[i + 1]
            concentration = None
            if pd.notna(raw_value) and str(raw_value).strip() != "":
                try:
                    concentration = float(raw_value)
                except ValueError:
                    concentration = None

            cursor.execute(
                "INSERT INTO results (sample_id, analyte_id, concentration) VALUES (?, ?, ?)",
                (sample_id, analyte_id_map[analyte_name], concentration),
            )

        imported_count += 1

    conn.commit()
    conn.close()
    print(f"Imported {imported_count} samples from {meta['source_filename']} (run_date={meta['run_date']}, panel={meta['panel']})")


# ==============================================================================
# CSV IMPORTER — NEW FORMAT
# ==============================================================================


def import_csv_new(csv_path: str, db_path=None):
    """Import new-format CSV with full metadata columns."""
    conn = get_connection(db_path)
    cursor = conn.cursor()

    meta = parse_filename_new(csv_path)

    cursor.execute("SELECT run_id FROM runs WHERE source_filename = ?", (meta["source_filename"],))
    if cursor.fetchone():
        print(f"Already imported: {meta['source_filename']}")
        conn.close()
        return

    # Read CSV: row 1 is the top header (Sample,,,...,analyte Results,...)
    # row 2 is sub-header (Name, Data File, Data Path, Type, Level, Acq. Date-Time, Sample Group, Pos., , Calc. Conc., ...)
    # Data starts at row 3
    df_raw = pd.read_csv(csv_path, header=None)

    # Row 0 has the top-level headers: find analyte columns
    top_headers = df_raw.iloc[0].tolist()

    # Find the first analyte Results column index
    analyte_start_idx = None
    analyte_columns = []
    for i, h in enumerate(top_headers):
        h_str = str(h).strip() if pd.notna(h) else ""
        if "Results" in h_str:
            if analyte_start_idx is None:
                analyte_start_idx = i
            analyte_columns.append(h_str.replace(" Results", "").strip())

    if analyte_start_idx is None:
        raise ValueError("Could not find analyte Results columns in CSV header")

    # Row 1 has sub-headers — identify metadata column positions
    sub_headers = df_raw.iloc[1].tolist()
    # Metadata columns are before analyte_start_idx
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

    # Extract data_path from first data row for the run record
    first_data_path = None
    if col_data_path is not None and len(df_raw) > 2:
        first_data_path = str(df_raw.iloc[2, col_data_path]).strip() if pd.notna(df_raw.iloc[2, col_data_path]) else None

    cursor.execute(
        "INSERT INTO runs (run_date, panel, source_filename, method_name, data_path) VALUES (?, ?, ?, ?, ?)",
        (meta["run_date"], meta["panel"], meta["source_filename"], meta["method_name"], first_data_path),
    )
    run_id = cursor.lastrowid

    # Lookup analyte IDs
    analyte_id_map = {}
    for name in analyte_columns:
        cursor.execute("SELECT analyte_id FROM analytes WHERE name = ?", (name,))
        row = cursor.fetchone()
        if row:
            analyte_id_map[name] = row[0]
        else:
            raise ValueError(f"Unknown analyte: {name}. Run 'init' first.")

    cursor.execute("SELECT type_code, type_id FROM sample_types")
    type_map = dict(cursor.fetchall())

    # Process data rows (skip header rows 0 and 1)
    imported_count = 0
    for row_idx in range(2, len(df_raw)):
        row = df_raw.iloc[row_idx]

        # Get data filename
        data_filename = None
        if col_data_file is not None:
            val = row.iloc[col_data_file]
            data_filename = str(val).strip() if pd.notna(val) else None

        if not data_filename or data_filename == "nan":
            continue

        # Get metadata
        sample_name = str(row.iloc[col_name]).strip() if col_name is not None and pd.notna(row.iloc[col_name]) else None
        instrument_type = str(row.iloc[col_type]).strip() if col_type is not None and pd.notna(row.iloc[col_type]) else None
        level = str(row.iloc[col_level]).strip() if col_level is not None and pd.notna(row.iloc[col_level]) else None
        acq_datetime = str(row.iloc[col_acq_datetime]).strip() if col_acq_datetime is not None and pd.notna(row.iloc[col_acq_datetime]) else None
        sample_group = str(row.iloc[col_sample_group]).strip() if col_sample_group is not None and pd.notna(row.iloc[col_sample_group]) else None
        position = str(row.iloc[col_pos]).strip() if col_pos is not None and pd.notna(row.iloc[col_pos]) else None

        # Classify sample using instrument type when available
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
            (
                run_id, data_filename, sample_name, sample_type_id, instrument_type,
                acq_datetime, position, sample_group,
                info.collection_date, info.patient_sequence,
                info.calibrator_level, info.qc_level, info.qc_replicate,
                info.eqa_scheme, info.eqa_year, info.eqa_round,
                info.eqa_sample_code, info.eqa_replicate,
            ),
        )
        sample_id = cursor.lastrowid

        # Insert concentration results
        for i, analyte_name in enumerate(analyte_columns):
            col_idx = analyte_start_idx + i
            raw_value = row.iloc[col_idx] if col_idx < len(row) else None
            concentration = None
            if pd.notna(raw_value) and str(raw_value).strip() != "":
                try:
                    concentration = float(raw_value)
                except ValueError:
                    concentration = None

            cursor.execute(
                "INSERT INTO results (sample_id, analyte_id, concentration) VALUES (?, ?, ?)",
                (sample_id, analyte_id_map[analyte_name], concentration),
            )

        imported_count += 1

    conn.commit()
    conn.close()
    print(f"Imported {imported_count} samples from {meta['source_filename']} (method={meta['method_name']}, run_date={meta['run_date']}, panel={meta['panel']})")


# ==============================================================================
# UNIFIED IMPORTER
# ==============================================================================


def import_csv(csv_path: str, db_path=None):
    """Auto-detect format and import CSV."""
    fmt = detect_format(csv_path)
    if fmt == "old":
        import_csv_old(csv_path, db_path)
    else:
        import_csv_new(csv_path, db_path)


# ==============================================================================
# QUERIES
# ==============================================================================


def query_qc_trend(analyte_name: str, qc_level: str, db_path=None):
    """Get QC results for Levey-Jennings charting."""
    conn = get_connection(db_path)
    query = """
        SELECT
            r.run_date,
            r.method_name,
            s.data_filename,
            s.qc_replicate,
            s.acquisition_datetime,
            res.concentration,
            qt.target_mean,
            qt.target_sd,
            CASE
                WHEN qt.target_sd > 0 THEN
                    (res.concentration - qt.target_mean) / qt.target_sd
                ELSE NULL
            END as z_score
        FROM results res
        JOIN samples s ON res.sample_id = s.sample_id
        JOIN runs r ON s.run_id = r.run_id
        JOIN analytes a ON res.analyte_id = a.analyte_id
        JOIN sample_types st ON s.sample_type_id = st.type_id
        LEFT JOIN qc_targets qt ON (
            qt.analyte_id = a.analyte_id
            AND qt.qc_level = s.qc_level
            AND qt.effective_from <= r.run_date
            AND (qt.effective_to IS NULL OR qt.effective_to >= r.run_date)
        )
        WHERE a.name = ?
          AND st.type_code = 'qc'
          AND s.qc_level = ?
          AND res.concentration IS NOT NULL
        ORDER BY r.run_date, s.qc_replicate
    """
    df = pd.read_sql_query(query, conn, params=[analyte_name, qc_level])
    conn.close()
    return df


def query_eqa_results(scheme: str, year: int, round_num: int, db_path=None):
    """Get all EQA results for a given scheme/year/round."""
    conn = get_connection(db_path)
    query = """
        SELECT
            a.name as analyte,
            s.eqa_sample_code,
            s.eqa_replicate,
            s.sample_group,
            s.acquisition_datetime,
            res.concentration,
            et.consensus_mean,
            et.consensus_sd,
            CASE
                WHEN et.consensus_sd > 0 THEN
                    (res.concentration - et.consensus_mean) / et.consensus_sd
                ELSE NULL
            END as z_score
        FROM results res
        JOIN samples s ON res.sample_id = s.sample_id
        JOIN runs r ON s.run_id = r.run_id
        JOIN analytes a ON res.analyte_id = a.analyte_id
        LEFT JOIN eqa_targets et ON (
            et.analyte_id = a.analyte_id
            AND et.scheme = s.eqa_scheme
            AND et.year = s.eqa_year
            AND et.round = s.eqa_round
            AND et.sample_code = s.eqa_sample_code
        )
        WHERE s.eqa_scheme = ?
          AND s.eqa_year = ?
          AND s.eqa_round = ?
        ORDER BY a.display_order, s.eqa_sample_code, s.eqa_replicate
    """
    df = pd.read_sql_query(query, conn, params=[scheme, year, round_num])
    conn.close()
    return df


def query_run_summary(db_path=None):
    """Get summary of all runs."""
    conn = get_connection(db_path)
    query = """
        SELECT
            r.run_id, r.run_date, r.panel, r.method_name, r.source_filename, r.imported_at,
            COUNT(DISTINCT s.sample_id) as sample_count
        FROM runs r
        LEFT JOIN samples s ON r.run_id = s.run_id
        GROUP BY r.run_id ORDER BY r.run_date DESC
    """
    df = pd.read_sql_query(query, conn)
    conn.close()
    return df


# ==============================================================================
# QC EVALUATION (Westgard Rules)
# ==============================================================================


def evaluate_westgard(values: list, mean: float, sd: float) -> list:
    """Evaluate Westgard multi-rules on a list of QC values."""
    violations = []
    for i, val in enumerate(values):
        z = (val - mean) / sd if sd > 0 else 0
        point_violations = []

        if abs(z) > 3.0:
            point_violations.append("1-3s")
        elif abs(z) > 2.0:
            point_violations.append("1-2s_warning")

        if i >= 1:
            z_prev = (values[i - 1] - mean) / sd if sd > 0 else 0
            if (z > 2.0 and z_prev > 2.0) or (z < -2.0 and z_prev < -2.0):
                point_violations.append("2-2s")
            if abs(z - z_prev) > 4.0:
                point_violations.append("R-4s")

        if i >= 3:
            recent = [(values[j] - mean) / sd for j in range(i - 3, i + 1)] if sd > 0 else []
            if recent and (all(v > 1.0 for v in recent) or all(v < -1.0 for v in recent)):
                point_violations.append("4-1s")

        if i >= 9:
            recent = [(values[j] - mean) / sd for j in range(i - 9, i + 1)] if sd > 0 else []
            if recent and (all(v > 0 for v in recent) or all(v < 0 for v in recent)):
                point_violations.append("10x")

        if point_violations:
            violations.append({"index": i, "value": val, "z_score": z, "rules": point_violations})

    return violations


# ==============================================================================
# CLI
# ==============================================================================


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return

    command = sys.argv[1]

    if command == "init":
        init_db()

    elif command == "import":
        if len(sys.argv) < 3:
            print("Usage: python steroid_panel_lims.py import <csv_path>")
            return
        csv_path = sys.argv[2]
        if not Path(csv_path).exists():
            print(f"File not found: {csv_path}")
            return
        import_csv(csv_path)

    elif command == "qc":
        if len(sys.argv) < 4:
            print("Usage: python steroid_panel_lims.py qc <analyte_name> <Low|High>")
            return
        analyte = sys.argv[2]
        level = sys.argv[3]
        df = query_qc_trend(analyte, level)
        if df.empty:
            print(f"No QC data found for {analyte} ({level})")
        else:
            print(df.to_string(index=False))

    elif command == "eqa":
        if len(sys.argv) < 5:
            print("Usage: python steroid_panel_lims.py eqa <scheme> <year> <round>")
            return
        scheme = sys.argv[2]
        year = int(sys.argv[3])
        round_num = int(sys.argv[4])
        df = query_eqa_results(scheme, year, round_num)
        if df.empty:
            print(f"No EQA data found for {scheme} {year} round {round_num}")
        else:
            print(df.to_string(index=False))

    elif command == "runs":
        df = query_run_summary()
        if df.empty:
            print("No runs imported yet.")
        else:
            print(df.to_string(index=False))

    else:
        print(f"Unknown command: {command}")
        print(__doc__)


if __name__ == "__main__":
    main()
