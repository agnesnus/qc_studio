# QC Studio

QC Studio is a Streamlit application for:
- QC data import (CSV and Excel)
- run storage in SQLite
- dashboard charts with QC target mean and SD tracking
- export and report views

---

## Quick Start

1. Install dependencies

```bash
pip install -r requirements.txt
```

2. Run the app

```bash
streamlit run qc_unified_app.py
```

---

## Reuse This Codebase For A New Project

Use this repository as a base template when you want to build another Streamlit website with a similar workflow.

### 1. Duplicate the project

Option A: Use this repository as a template on GitHub.

Option B: Clone locally and copy the folder:

```bash
git clone https://github.com/agnesnus/qc_studio.git
cp -R qc_studio your_new_project
cd your_new_project
```

### 2. Rename the app entry file (optional)

If you want a new app name:

```bash
mv qc_unified_app.py app.py
```

Then run with:

```bash
streamlit run app.py
```

### 3. Update project identity

Edit these first:
- `README.md` title and description
- app title/header in `qc_unified_app.py`
- module names or menu labels (`Dashboard`, `Database`, `Export`, `Report`)

### 4. Define your new data model

For a new domain/project, update:
- `ANALYTES` list (or equivalent entities)
- `SAMPLE_TYPES` if needed
- SQL schema (`SCHEMA_SQL`) tables/fields
- importer mappings and column detection logic

### 5. Update import pipelines

Modify these functions for your new input format:
- `import_csv_old`
- `import_csv_new`
- `import_excel_qc_file`
- helper parsers (filename parsing, column detection, normalization)

Tip: keep one import route stable first (for example CSV), then add complex Excel parsing.

### 6. Update chart and report logic

Adjust these parts:
- `make_qc_chart` for your new metrics
- `generate_final_report` for your summary columns
- `create_value_pictogram` if you want a different visual indicator

### 7. Keep target versioning behavior (recommended)

If your new project also has changing reference ranges over time, keep:
- `qc_targets` table design with `effective_from` and `effective_to`
- `get_qc_target` and `get_per_date_targets`
- step-line rendering in charts

This enables date-aware reference values without rewriting historical records.

### 8. Add target bulk upload for your new project

The app already supports:
- manual target entry
- bulk target upload from CSV/Excel

You can reuse `import_qc_targets_file` directly and just rename expected columns if needed.

### 9. Validate after each change

```bash
python3 -m py_compile qc_unified_app.py
streamlit run qc_unified_app.py
```

---

## Suggested Structure For New Projects

- `app.py` or `qc_unified_app.py`: main Streamlit app
- `requirements.txt`: dependencies
- `.streamlit/config.toml`: theme/runtime settings
- `data/` or test folders for sample uploads
- SQLite DB file generated at runtime

---

## Deployment Notes

For Streamlit Community Cloud:
- set main file to `qc_unified_app.py` (or your renamed entry file)
- ensure `requirements.txt` includes all used packages
- commit and push changes to GitHub

---

## Adaptation Checklist

- [ ] New project name and branding updated
- [ ] Data schema updated
- [ ] Import logic updated for new file formats
- [ ] Dashboard and report metrics updated
- [ ] Target/reference range logic verified
- [ ] README updated with domain-specific instructions
- [ ] Smoke test done locally

---

## License

See `LICENSE`.
