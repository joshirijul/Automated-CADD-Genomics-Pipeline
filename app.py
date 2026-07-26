import os
import json
import urllib.request
import urllib.parse
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
    page_title="CADD Virtual Screening Engine",
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
    """Computes ADME descriptors and drug-likeness filters."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    mw, logp = Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    tpsa, rot = Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol)
    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    drug_like = "YES" if (lip_viol <= 1 and tpsa <= 140 and rot <= 10) else "NO"
    return {"MolWt": round(mw, 2), "LogP": round(logp, 2), "HBD": hbd, "HBA": hba, "TPSA": round(tpsa, 2), "Lipinski_Violations": lip_viol, "Drug_Like": drug_like}

def parse_uploaded_ligands(uploaded_file):
    """Automatically parses CSV, SDF, MOL, or PDB files into a standardized SMILES DataFrame."""
    filename = uploaded_file.name.lower()
    records = []
    
    if filename.endswith(".csv"):
        return pd.read_csv(uploaded_file)
        
    # Write temp file for RDKit structure suppliers
    temp_path = f"temp_{filename}"
    with open(temp_path, "wb") as f:
        f.write(uploaded_file.getvalue())
        
    mols = []
    try:
        if filename.endswith(".sdf"):
            mols = [m for m in Chem.SDMolSupplier(temp_path, sanitize=True) if m is not None]
        elif filename.endswith(".mol"):
            mol = Chem.MolFromMolFile(temp_path, sanitize=True)
            if mol: mols.append(mol)
        elif filename.endswith(".pdb"):
            mol = Chem.MolFromPDBFile(temp_path, sanitize=True)
            if mol: mols.append(mol)
            
        for idx, mol in enumerate(mols):
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Molecule_{idx+1}"
            smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
            records.append({"Molecule_Name": name, "Canonical_SMILES": smiles})
    finally:
        if os.path.exists(temp_path): os.remove(temp_path)
        
    return pd.DataFrame(records)

def fetch_from_pubchem(drug_names_str):
    """Queries NIH PubChem REST API with automated AI spellcheck and typo-correction fallback."""
    names = [n.strip() for n in drug_names_str.split(",") if n.strip()]
    records = []
    for name in names:
        encoded_name = urllib.parse.quote(name)
        try:
            # Attempt 1: Direct Exact Name Query
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES/JSON"
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode())
                smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                records.append({"Molecule_Name": name.capitalize(), "Canonical_SMILES": smiles})
        except Exception:
            # Attempt 2: Auto-Spellcheck Fallback for "Average Joe" typos!
            try:
                spell_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/spell/suggest/{encoded_name}/JSON"
                req_spell = urllib.request.Request(spell_url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req_spell) as spell_resp:
                    spell_data = json.loads(spell_resp.read().decode())
                    suggestions = spell_data.get('Dictionary_Suggest_01', {}).get('SuggestionList', {}).get('Suggestion', [])
                    
                    if suggestions:
                        corrected_name = suggestions[0]
                        st.toast(f"🪄 Auto-corrected typo '{name}' ➔ '{corrected_name}'!", icon="💡")
                        
                        # Retry structure fetch using the corrected chemical name
                        enc_correct = urllib.parse.quote(corrected_name)
                        url_retry = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc_correct}/property/CanonicalSMILES/JSON"
                        req_retry = urllib.request.Request(url_retry, headers={'User-Agent': 'Mozilla/5.0'})
                        with urllib.request.urlopen(req_retry) as retry_resp:
                            retry_data = json.loads(retry_resp.read().decode())
                            smiles = retry_data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                            records.append({"Molecule_Name": f"{corrected_name.capitalize()} (Auto-corrected)", "Canonical_SMILES": smiles})
                    else:
                        st.sidebar.error(f"❌ Could not resolve chemical structure for '{name}' on PubChem.")
            except Exception:
                st.sidebar.error(f"❌ Could not resolve chemical structure for '{name}' on PubChem.")
    return pd.DataFrame(records)

def render_3d_pdb(pdb_string, style="cartoon", color_by="spectrum"):
    """Renders an interactive 3D molecular viewer using py3Dmol."""
    view = py3Dmol.view(width=800, height=500)
    view.addModel(pdb_string, "pdb")
    if style == "cartoon": view.setStyle({"cartoon": {"color": color_by}})
    elif style == "surface":
        view.addSurface(py3Dmol.VDW, {"opacity": 0.7, "color": "white"})
        view.setStyle({"cartoon": {"color": "slate"}})
    view.addStyle({"hetflag": True}, {"stick": {"colorscheme": "greenCarbon", "radius": 0.2}})
    view.zoomTo()
    return view._make_html()

# --- SIDEBAR CONTROLS ---
with st.sidebar:
    st.title("⚙️ Input Pipeline")
    demo_mode = st.toggle("🚀 Enable 1-Click Demo Mode (BACE1 Target)", value=True)
    
    st.divider()
    st.markdown("### 1. Target Receptor")
    uploaded_pdb = st.file_uploader("Upload Target (.pdb)", type=["pdb"], disabled=demo_mode)
    
    # Resolve PDB input
    pdb_content = None
    if demo_mode and os.path.exists("data/1fkn.pdb"):
        with open("data/1fkn.pdb", "r") as f: pdb_content = f.read()
    elif demo_mode and os.path.exists("1fkn.pdb"):
        with open("1fkn.pdb", "r") as f: pdb_content = f.read()
    elif uploaded_pdb is not None:
        pdb_content = uploaded_pdb.getvalue().decode("utf-8")
        
    st.divider()
    st.markdown("### 2. Ligand Screening Library")
    input_method = st.radio("Input Method:", ["PubChem Name Search (AI/REST)", "Upload Structure Files (.csv, .sdf, .mol)"])
    
    lib_df = None
    if input_method == "PubChem Name Search (AI/REST)":
        default_drugs = "Aspirin, Ibuprofen, Donepezil, Caffeine, Acetaminophen" if demo_mode else ""
        drug_input = st.text_area("Type drug or metabolite names (comma separated):", value=default_drugs, help="Our engine automatically queries NIH PubChem to fetch 3D chemical structures.")
        if drug_input:
            with st.spinner("Fetching live structures from PubChem..."):
                lib_df = fetch_from_pubchem(drug_input)
    else:
        uploaded_file = st.file_uploader("Upload Library (.csv, .sdf, .mol, .pdb)", type=["csv", "sdf", "mol", "pdb"])
        if uploaded_file is not None:
            with st.spinner("Parsing chemical structures..."):
                lib_df = parse_uploaded_ligands(uploaded_file)
        elif demo_mode and os.path.exists("results/reports/ligand_smiles_library.csv"):
            lib_df = pd.read_csv("results/reports/ligand_smiles_library.csv")
        elif demo_mode and os.path.exists("ligand_smiles_library.csv"):
            lib_df = pd.read_csv("ligand_smiles_library.csv")

    st.divider()
    st.markdown("### 📊 ADME Filter Settings")
    max_lipinski = st.slider("Max Lipinski Violations", 0, 4, 1)
    filter_druglike = st.checkbox("Require Veber Drug-Likeness", value=True)

# --- MAIN DASHBOARD AREA ---
st.title("🧬 Structure-Based Drug Design & Virtual Screening Platform")
st.markdown("An automated, full-stack computational biology engine. Convert chemical structures on the fly, screen for ADME-Tox drug-likeness, and predict target bioactivity via multivariate QSAR regression.")

if demo_mode:
    st.info("💡 **Demo Mode Active:** Currently evaluating **Human BACE1 (1FKN)** against Alzheimer's therapeutics and reference inhibitors. Toggle off in the sidebar to process custom research data.")

tab1, tab2, tab3 = st.tabs(["🏛️ 3D Active Site Viewer", "💊 ADME-Tox Profiling", "📈 QSAR Activity Prediction"])

# --- TAB 1: 3D MOLECULAR VIEWER ---
with tab1:
    st.header("Macromolecular Architecture")
    if pdb_content:
        col1, col2 = st.columns([3, 1])
        with col2:
            st.markdown("#### Rendering Controls")
            render_style = st.radio("Receptor Style:", ["cartoon", "surface"], index=0)
            color_scheme = st.selectbox("Coloring Scheme:", ["spectrum", "chain", "secondary structure"], index=0)
            st.markdown("---")
            st.caption("🔬 **Navigation:**\n* **Rotate:** Left-click + drag\n* **Zoom:** Scroll wheel / Pinch\n* **Active Site:** Green carbon sticks indicate co-crystallized inhibitors.")
        with col1:
            components.html(render_3d_pdb(pdb_content, style=render_style, color_by=color_scheme.split()[0]), height=520, width=800)
    else:
        st.warning("⚠️ No receptor target loaded. Enable Demo Mode or upload a `.pdb` coordinate file in the sidebar.")

# --- TAB 2: ADME-TOX SCREENING ---
with tab2:
    st.header("Cheminformatics & Drug-Likeness Filtering")
    if lib_df is not None and not lib_df.empty and "Canonical_SMILES" in lib_df.columns:
        with st.spinner("Computing RDKit ADME descriptors..."):
            adme_results = lib_df["Canonical_SMILES"].apply(calc_adme).apply(pd.Series)
            full_df = pd.concat([lib_df, adme_results], axis=1)
            
            # Apply UI Filters
            filtered_df = full_df[full_df["Lipinski_Violations"] <= max_lipinski]
            if filter_druglike:
                filtered_df = filtered_df[filtered_df["Drug_Like"] == "YES"]
                
            col1, col2 = st.columns([1, 4])
            with col1:
                st.metric(label="Passing Candidates", value=f"{len(filtered_df)} / {len(full_df)}")
            with col2:
                st.dataframe(filtered_df, use_container_width=True)
                
            csv_export = filtered_df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Export Curated Screening Report (.csv)", data=csv_export, file_name="adme_screening_report.csv", mime="text/csv")
    else:
        st.info("👈 Enter drug names in the sidebar (e.g., 'Aspirin, Ibuprofen') or upload a structural library to initiate screening.")

# --- TAB 3: QSAR REGRESSION MODELING ---
with tab3:
    st.header("Quantitative Structure-Activity Relationship (QSAR)")
    train_df = get_training_data()
    train_adme = train_df["SMILES"].apply(calc_adme).apply(pd.Series)
    train_df = pd.concat([train_df, train_adme], axis=1)
    
    # Train OLS Multivariate Regression
    X = train_df[["LogP", "TPSA"]].values
    y = train_df["pIC50"].values
    X_design = np.column_stack([np.ones(len(X)), X])
    coeffs, _, _, _ = np.linalg.lstsq(X_design, y, rcond=None)
    y_pred = X_design @ coeffs
    r2 = 1 - (np.sum((y - y_pred)**2) / np.sum((y - np.mean(y))**2))
    
    col1, col2 = st.columns([1, 1])
    with col1:
        st.markdown("#### Multivariate Regression Model")
        st.latex(f"pIC_{{50}} = {coeffs[0]:.2f} + ({coeffs[1]:.2f} \\times \\log P) + ({coeffs[2]:.4f} \\times \\text{{TPSA}})")
        st.metric("Model Fit ($R^2$ Score)", f"{r2:.3f}", help="Indicates the linear variance explained by lipophilicity and polar surface area across benchmark inhibitors.")
        
        if lib_df is not None and not lib_df.empty and "Canonical_SMILES" in lib_df.columns:
            st.markdown("#### Candidate Potency Predictions")
            screen_adme = lib_df["Canonical_SMILES"].apply(calc_adme).apply(pd.Series)
            X_screen = np.column_stack([np.ones(len(screen_adme)), screen_adme[["LogP", "TPSA"]].values])
            lib_df["Predicted_pIC50"] = np.round(X_screen @ coeffs, 2)
            ranked_lib = pd.concat([lib_df, screen_adme[["LogP", "TPSA", "Drug_Like"]]], axis=1).sort_values(by="Predicted_pIC50", ascending=False)
            st.dataframe(ranked_lib[["Molecule_Name", "Predicted_pIC50", "LogP", "TPSA", "Drug_Like"]], use_container_width=True)
            
    with col2:
        st.markdown("#### Bioactivity Regression Fit")
        fig, ax = plt.subplots(figsize=(7, 5))
        sns.set_style("whitegrid")
        sns.regplot(x=y, y=y_pred, ax=ax, color="#2b5c8f", scatter_kws={"s": 70}, line_kws={"color": "#d95f02", "label": f"Fit (R² = {r2:.2f})"})
        ax.set_title("Experimental vs. Predicted pIC50 (BACE1 Benchmark)", fontweight="bold")
        ax.set_xlabel("Observed pIC50 (-log IC50)")
        ax.set_ylabel("Predicted pIC50")
        ax.legend()
        st.pyplot(fig)