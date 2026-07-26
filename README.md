# 🧬 Automated CADD Pre-Processing & Virtual Screening Pipeline

An enterprise-grade, hybrid computational biology platform designed for structure-based (SBDD) and ligand-based (LBDD) drug discovery. 

This repository bridges the gap between **heavy command-line automation** (for CPU-intensive 3D docking and thermodynamic simulations) and an **interactive web dashboard** (for rapid macromolecular inspection, multi-format chemical standardization, and ADME-Tox library curation).

---

## 🏛️ System Architecture: CLI vs. Web Dashboard

Running 3D grid-search molecular docking (AutoDock Vina) and solvated molecular dynamics minimization (GROMACS) requires intensive, multi-hour CPU compute that cannot run reliably inside a standard web browser. Therefore, this platform is divided into two specialized tiers:

1. **Interactive Web Dashboard (`app.py`):** The pre-processing gateway. Used by bench scientists and researchers to visually validate 3D protein crystal structures, standardize raw ligand files, and filter out toxic candidates via cheminformatics rules *before* expending supercomputer CPU hours.
2. **Command-Line Backend (`scripts/`):** The heavy-lifting automation engine. Executes high-throughput sequence alignment, active-site topography mapping, automated docking grid generation, and GROMACS thermodynamic system solvation.

---

## 💻 Full-Stack Interactive Workbench (`app.py` v3.0)

A reactive **Streamlit** web application featuring 3D macromolecular rendering via **py3Dmol** and 2D cheminformatics curation via **RDKit**.

### Key Features:
* **🏛️ 3D Macromolecule Inspector (SBDD):** Upload any receptor `.pdb` coordinate file for interactive 3D rotation, zooming, surface rendering, and active-site verification directly in the browser.
* **💊 Automated Library Curation & ADME-Tox (LBDD):** Standardize raw multi-format chemical libraries (`.sdf`, `.mol`, `.csv`, `.pdb`) on the fly. Computes real-time RDKit descriptors (Molecular Weight, $\log P$, TPSA, Rotatable Bonds) to screen candidates against **Lipinski's Rule of 5** and **Veber's Drug-Likeness Criteria**.
* **🪄 3-Tier Fault-Tolerant PubChem Fetcher:**
  1. *Instant Reference Cache:* Bypasses cloud datacenter IP rate-limiting for benchmark therapeutics.
  2. *Live REST API:* Queries NIH PubChem using NCBI-compliant academic research headers.
  3. *AI Spellcheck Fallback:* Automatically detects user typos in chemical names (e.g., `"Caffiene"` $\rightarrow$ `"Caffeine"`), queries the PubChem spelling suggestion engine, and resolves the corrected 3D structure seamlessly.
* **📥 One-Click Export:** Download curated, drug-like screening tables as clean CSVs ready for immediate command-line docking ingestion.

### Launch the Dashboard Locally:
```bash
conda activate cadd_env
streamlit run app.py