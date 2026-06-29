# QC Studio

This project is a ready-made website template built with Streamlit.

You can use it as a base and create a new website for a different project, even if you have never coded before.

---

## What this app already does

- Lets users upload files
- Stores data in a local database
- Shows charts and summaries
- Exports results

You can reuse this structure and replace the content with your own project topic.

---

## Before you start

You only need these things:

1. A GitHub account
2. A computer with internet
3. Python installed (version 3.10 or newer is recommended)
4. A text editor (VS Code is easiest)

---

## Super Simple Plan

1. Copy this project
2. Change the app title and text
3. Replace the data fields with your project fields
4. Test locally
5. Publish online

---

## Step-by-step for complete beginners

### Step 1: Copy this project to your own GitHub

- Open this repository on GitHub
- Click Use this template (or create a new repository and upload these files)
- Give your new repository a name

Now you have your own copy.

### Step 2: Download your copy to your computer

If you know Git, clone it.

If you do not know Git:
- Click Code on GitHub
- Click Download ZIP
- Unzip it on your computer

### Step 3: Open the folder in VS Code

- Open VS Code
- Choose File > Open Folder
- Select your project folder

### Step 4: Install required packages

Open the terminal in VS Code and run:

```bash
pip install -r requirements.txt
```

### Step 5: Start the website locally

Run:

```bash
streamlit run qc_unified_app.py
```

A browser page will open. This is your website running on your computer.

### Step 6: Change the website text first

In file qc_unified_app.py, change:
- Main title
- Section names
- Help text

This is the easiest first change and helps you learn the file layout.

### Step 7: Replace project fields with your own

This app currently uses QC terms like analyte, HQC, and LQC.

For your project, replace these with your own terms.

Example:
- analyte -> product
- qc_level -> category
- concentration -> score

Update in these areas:
- Data lists near the top of the file
- Table schema (database section)
- Upload/import logic
- Charts and report labels

### Step 8: Keep what you do not need to rebuild

You can keep these parts and just rename fields:
- Upload flow
- Database saving
- Chart layout
- Report table
- CSV export

This saves a lot of time.

### Step 9: Test after each small change

After edits, run:

```bash
python3 -m py_compile qc_unified_app.py
streamlit run qc_unified_app.py
```

If the page opens and works, continue.

---

## Publish your website online (easy method)

Use Streamlit Community Cloud.

1. Push your code to GitHub
2. Go to share.streamlit.io
3. Sign in with GitHub
4. Click New app
5. Choose your repository
6. Set main file to qc_unified_app.py
7. Click Deploy

Your website is now online with a public link.

---

## If you want to build a different project quickly

Use this checklist:

- [ ] Rename app title and sections
- [ ] Replace data fields with your project fields
- [ ] Update upload file format rules
- [ ] Update chart labels and report columns
- [ ] Test locally
- [ ] Deploy to Streamlit Cloud

---

## Common beginner tips

- Change one thing at a time
- Test after every change
- Keep a backup copy before big edits
- If something breaks, undo only your last small change

---

## Need help

If you get stuck, open an issue in the repository and include:
- What you changed
- The error message
- A screenshot

That makes troubleshooting much faster.

---

## License

See LICENSE.
