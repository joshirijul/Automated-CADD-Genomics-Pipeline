import os
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import py3Dmol

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="CADD Virtual Screening Dashboard",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- BACKEND LOGIC & DATA LOADERS ---
@st.cache_data
def get_training_data():
    """Returns the benchmark BACE1 training set for QSAR modeling."""
    data = [
        {"Name": "BACE_Inh_1", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)C#N)F", "pIC50": 7.82},
        {"Name": "BACE_Inh_2", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)C(=O)N)c4cncnc4", "pIC50": 8.15},
        {"Name": "BACE_Inh_3", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)OC)c4cccnc4", "pIC50": 7.40},
        {"Name": "BACE_Inh_4", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)C(=O)NC)c4cncnc4", "pIC50": 8.30},
        {"Name": "BACE_Inh_5", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)Cl)c4cncnc4", "pIC50": 7.10},
        {"Name": "BACE_Inh_6", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)F)F", "pIC50": 6.50},
        {"Name": "BACE_Inh_7", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)N(C)C)F", "pIC50": 6.85},
        {"Name": "BACE_Inh_8", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)S(=O)(=O)C)c4cncnc4", "pIC50": 7.95},
        {"Name": "BACE_Inh_9", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)C#N)c4cncnc4", "pIC50": 8.50},
        {"Name": "BACE_Inh_10", "SMILES": "CC1(C)c2cc(cc(c2N=C1N)c3cncc(c3)C)F", "pIC50": 6.20}
    ]
    return pd.DataFrame(data)

def calc_adme(smiles):
    """Computes ADME descriptors for a given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return None
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = Descriptors.TPSA(mol)
    rot = Descriptors.NumRotatableBonds(mol)
    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    drug_like = "YES" if (lip_viol <= 1 and tpsa <= 140 and rot <= 10) else "NO"
    return {"MolWt": round(mw, 2), "LogP": round(logp, 2), "HBD": hbd, "HBA": hba, "TPSA": round(tpsa, 2), "Lipinski_Violations": lip_viol, "Drug_Like": drug_like}

def render_3d_pdb(pdb_string, style="cartoon", color_by="spectrum"):
    """Renders an interactive 3D molecular viewer using py3Dmol."""
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_string, "pdb")
    if style == "cartoon":
        view.setStyle({"cartoon": {"color": color_by}})
    elif style == "surface":
        view.addSurface(py3Dmol.VDW, {"opacity": 0.7, "color": "white"})
        view.setStyle({"cartoon": {"color": "slate"}})
    # Highlight non-protein heteroatoms (bound ligands) as sticks
    view.addStyle({"hetflag": True}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}})
    view.zoomTo()
    return view._make_html()

# --- SIDEBAR CONTROLS ---
with st.sidebar:
    st.title("⚙️ Pipeline Controls")
    st.markdown("Upload your experimental structures to run automated virtual screening.")
    
    # File Uploader for Target PDB
    uploaded_pdb = st.file_uploader("1. Upload Target Receptor (.pdb)", type=["pdb"])
    
    # Load default PDB if none uploaded
    pdb_content = None
    if uploaded_pdb is not None:
        pdb_content = uploaded_pdb.getvalue().decode("utf-8")
    elif os.path.exists("data/1fkn.pdb"):
        with open("data/1fkn.pdb", "r") as f:
            pdb_content = f.read()
    elif os.path.exists("1fkn.pdb"):
        with open("1fkn.pdb", "r") as f:
            pdb_content = f.read()
            
    # File Uploader for Ligand Library
    uploaded_csv = st.file_uploader("2. Upload Ligand SMILES Library (.csv)", type=["csv"])
    
    st.divider()
    st.markdown("### 📊 Screening Parameters")
    max_lipinski = st.slider("Max Lipinski Violations Allowed", 0, 4, 1)
    filter_druglike = st.checkbox("Only show Drug-Like candidates (Veber Rules)", value=True)

# --- MAIN DASHBOARD AREA ---
st.title("🧬 Structure-Based Drug Design & Virtual Screening Engine")
st.markdown("Automated bioinformatics pipeline integrating **3D macromolecular visualization**, **ADME-Tox cheminformatics**, and **OLS QSAR bioactivity prediction**.")

tab1, tab2, tab3 = st.tabs(["🏛️ 3D Active Site Viewer", "💊 ADME-Tox & Drug-Likeness", "📈 QSAR Bioactivity Prediction"])

# --- TAB 1: 3D MOLECULAR VIEWER ---
with tab1:
    st.header("Interactive Target Macromolecule")
    if pdb_content:
        col1, col2 = st.columns([3, 1])
        with col2:
            st.markdown("#### Display Styles")
            render_style = st.radio("Receptor Representation:", ["cartoon", "surface"], index=0)
            color_scheme = st.selectbox("Color Scheme:", ["spectrum", "chain", "secondary structure"], index=0)
            st.info("💡 **Interactive Guide:**\n* **Left-click + Drag:** Rotate 3D structure.\n* **Scroll Wheel:** Zoom in/out.\n* **Green Sticks:** Co-crystallized inhibitor.")
        with col1:
            html_viewer = render_3d_pdb(pdb_content, style=render_style, color_by=color_scheme.split()[0])
            components.html(html_viewer, height=520, width=800)
    else:
        st.warning("⚠️ No PDB file detected. Please upload a target `.pdb` file in the sidebar or place `1fkn.pdb` in your project folder.")

# --- TAB 2: ADME-TOX SCREENING ---
with tab2:
    st.header("Chemical Library Profiling")
    lib_df = None
    if uploaded_csv is not None:
        lib_df = pd.read_csv(uploaded_csv)
    elif os.path.exists("results/reports/ligand_smiles_library.csv"):
        lib_df = pd.read_csv("results/reports/ligand_smiles_library.csv")
    elif os.path.exists("ligand_smiles_library.csv"):
        lib_df = pd.read_csv("ligand_smiles_library.csv")
        
    if lib_df is not None and "Canonical_SMILES" in lib_df.columns:
        # Calculate ADME on the fly
        with st.spinner("Calculating RDKit molecular descriptors..."):
            adme_results = lib_df["Canonical_SMILES"].apply(calc_adme).apply(pd.Series)
            full_df = pd.concat([lib_df, adme_results], axis=1)
            
            # Apply dynamic sidebar filters
            filtered_df = full_df[full_df["Lipinski_Violations"] <= max_lipinski]
            if filter_druglike:
                filtered_df = filtered_df[filtered_df["Drug_Like"] == "YES"]
                
            st.metric(label="Passing Candidates", value=f"{len(filtered_df)} / {len(full_df)}")
            st.dataframe(filtered_df, use_container_width=True)
            
            # Export button
            csv_export = filtered_df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Download Filtered Screening Report (.csv)", data=csv_export, file_name="filtered_adme_report.csv", mime="text/csv")
    else:
        st.info("📁 Please upload a valid CSV library with a `Canonical_SMILES` column in the sidebar to run ADME profiling.")

# --- TAB 3: QSAR REGRESSION MODELING ---
with tab3:
    st.header("Quantitative Structure-Activity Relationship (QSAR)")
    train_df = get_training_data()
    train_adme = train_df["SMILES"].apply(calc_adme).apply(pd.Series)
    train_df = pd.concat([train_df, train_adme], axis=1)
    
    # Train multivariate OLS Regression
    X = train_df[["LogP", "TPSA"]].values
    y = train_df["pIC50"].values
    X_design = np.column_stack([np.ones(len(X)), X])
    coeffs, _, _, _ = np.linalg.lstsq(X_design, y, rcond=None)
    y_pred = X_design @ coeffs
    
    # R^2 calculation
    r2 = 1 - (np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2))
    
    col1, col2 = st.columns([1, 1])
    with col1:
        st.markdown("#### Mathematical Model")
        st.latex(f"pIC_{{50}} = {coeffs[0]:.2f} + ({coeffs[1]:.2f} \\times \\log P) + ({coeffs[2]:.4f} \\times \\text{{TPSA}})")
        st.metric("Model Convergence ($R^2$ Score)", f"{r2:.3f}")
        
        if lib_df is not None and "Canonical_SMILES" in lib_df.columns:
            st.markdown("#### Workspace Candidate Predictions")
            screen_adme = lib_df["Canonical_SMILES"].apply(calc_adme).apply(pd.Series)
            X_screen = np.column_stack([np.ones(len(screen_adme)), screen_adme[["LogP", "TPSA"]].values])
            lib_df["Predicted_pIC50"] = np.round(X_screen @ coeffs, 2)
            ranked_lib = pd.concat([lib_df, screen_adme[["LogP", "TPSA", "Drug_Like"]]], axis=1).sort_values(by="Predicted_pIC50", ascending=False)
            st.dataframe(ranked_lib[["Molecule_Name", "Predicted_pIC50", "LogP", "TPSA", "Drug_Like"]], use_container_width=True)
            
    with col2:
        st.markdown("#### Regression Regression Analysis")
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.set_style("whitegrid")
        sns.regplot(x=y, y=y_pred, ax=ax, color="#2b5c8f", scatter_kws={"s": 70}, line_kws={"color": "#d95f02", "label": f"Fit (R² = {r2:.2f})"})
        ax.set_title("Experimental vs. Predicted pIC50", fontweight="bold")
        ax.set_xlabel("Observed pIC50")
        ax.set_ylabel("Predicted pIC50")
        ax.legend()
        st.pyplot(fig)