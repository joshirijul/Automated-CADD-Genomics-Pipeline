# Automated CADD & Genomics Pipeline: Targeted Virtual Screening & Biomolecular Characterization

![Python](https://img.shields.io/badge/Python-3.10%2B-blue?logo=python&logoColor=white)
![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-green)
![BioPython](https://img.shields.io/badge/BioPython-1.80%2B-yellow)
![GROMACS](https://img.shields.io/badge/GROMACS-Molecular_Dynamics-orange)
![AutoDock Vina](https://img.shields.io/badge/AutoDock_Vina-Molecular_Docking-purple)

An automated, deterministic, end-to-end computational biology and computer-aided drug design (CADD) pipeline. This repository bridges small-molecule cheminformatics, structural validation, QSAR regression modeling, molecular docking, and thermodynamic system minimization into a modular software architecture.

---

## 📂 Repository Architecture

```text
├── data/                    # Target macromolecule coordinates (1FKN) & raw ligand libraries (.sdf/.mol)
├── scripts/                 # Python automation engines & headless PyMol rendering macros
├── results/
│   ├── alignments/          # Multiple sequence alignments (ClustalW/MUSCLE) & phylogenetic trees
│   ├── docking/             # AutoDock Vina grid box configurations & prepared .pdbqt structures
│   ├── md_simulation/       # Solvated periodic box coordinates, topologies, and minimization XVG trajectories
│   └── reports/             # Curated QSAR datasets, chemical screening CSVs, and publication graphics
├── environment.yml          # Self-contained Conda environment dependencies
└── README.md                # Pipeline documentation