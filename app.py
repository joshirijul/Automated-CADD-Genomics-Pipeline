import os
import json
import urllib.request
import urllib.parse
import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
import py3Dmol

# --- PAGE CONFIGURATION ---
st.set_page_config(
    page_title="CADD Pre-Processing & Curation Workbench",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# --- BACKEND LOGIC & DATA LOADERS ---
def calc_adme(smiles):
    """Computes ADME descriptors and drug-likeness filters for any chemical ligand."""
    mol = Chem.MolFromSmiles(smiles)
    if not mol: return None
    mw, logp = Descriptors.MolWt(mol), Descriptors.MolLogP(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    tpsa, rot = Descriptors.TPSA(mol), Descriptors.NumRotatableBonds(mol)
    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    drug_like = "YES" if (lip_viol <= 1 and tpsa <= 140 and rot <= 10) else "NO"
    return {"MolWt": round(mw, 2), "LogP": round(logp, 2), "HBD": hbd, "HBA": hba, "TPSA": round(tpsa, 2), "Lipinski_Violations": lip_viol, "Drug_Like": drug_like}

def parse_uploaded_ligands(uploaded_file):
    """Parses CSV, SDF, MOL, or PDB ligand structure files into a standardized SMILES DataFrame."""
    filename = uploaded_file.name.lower()
    records = []
    
    if filename.endswith(".csv"):
        return pd.read_csv(uploaded_file)
        
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
    """
    3-Tier Fault-Tolerant Chemical Structure Fetcher:
    1. Instant Reference Cache (Bypasses cloud IP rate-limiting for common therapeutics)
    2. Live PubChem REST API (NCBI academic research headers)
    3. AI Spellcheck & Typo-Correction Fallback
    """
    BENCHMARK_SMILES = {
        "aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "donepezil": "COC1=C(C=C2C(=C1)CC(C2=O)CC3CCN(CC3)CC4=CC=CC=C4)OC",
        "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "acetaminophen": "CC(=O)NC1=CC=C(O)C=C1",
        "paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
        "metformin": "CN(C)C(=N)NC(=N)N",
        "atorvastatin": "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        "omeprazole": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
        "imatinib": "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5"
    }

    names = [n.strip() for n in drug_names_str.split(",") if n.strip()]
    records = []
    headers = {'User-Agent': 'CADD-PreProcessing-Workbench/3.0 (Academic Research; mailto:lab@university.edu)'}
    
    for name in names:
        clean_name = name.lower()
        if clean_name in BENCHMARK_SMILES:
            records.append({"Molecule_Name": name.capitalize(), "Canonical_SMILES": BENCHMARK_SMILES[clean_name]})
            continue
            
        encoded_name = urllib.parse.quote(name)
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/property/CanonicalSMILES/JSON"
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req) as response:
                data = json.loads(response.read().decode())
                smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                records.append({"Molecule_Name": name.capitalize(), "Canonical_SMILES": smiles})
        except Exception:
            try:
                spell_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/spell/suggest/{encoded_name}/JSON"
                req_spell = urllib.request.Request(spell_url, headers=headers)
                with urllib.request.urlopen(req_spell) as spell_resp:
                    spell_data = json.loads(spell_resp.read().decode())
                    suggestions = spell_data.get('Dictionary_Suggest_01', {}).get('SuggestionList', {}).get('Suggestion', [])
                    
                    if suggestions:
                        corrected = suggestions[0]
                        if corrected.lower() in BENCHMARK_SMILES:
                            st.toast(f"🪄 Auto-corrected typo '{name}' ➔ '{corrected}'!", icon="💡")
                            records.append({"Molecule_Name": f"{corrected.capitalize()} (Auto-corrected)", "Canonical_SMILES": BENCHMARK_SMILES[corrected.lower()]})
                        else:
                            enc_correct = urllib.parse.quote(corrected)
                            url_retry = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{enc_correct}/property/CanonicalSMILES/JSON"
                            req_retry = urllib.request.Request(url_retry, headers=headers)
                            with urllib.request.urlopen(req_retry) as retry_resp:
                                retry_data = json.loads(retry_resp.read().decode())
                                smiles = retry_data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                                st.toast(f"🪄 Auto-corrected typo '{name}' ➔ '{corrected}'!", icon="💡")
                                records.append({"Molecule_Name": f"{corrected.capitalize()} (Auto-corrected)", "Canonical_SMILES": smiles})
                    else:
                        st.sidebar.error(f"❌ Could not resolve '{name}' on PubChem (Compound not found).")
            except Exception:
                st.sidebar.error(f"❌ Could not resolve '{name}'. Note: External API rate limits may apply on cloud servers.")
                
    return pd.DataFrame(records)

def render_3d_pdb(pdb_string, style="cartoon", color_by="spectrum"):
    """Renders an interactive 3D macromolecular viewer using py3Dmol."""
    view = py3Dmol.view(width=800, height=520)
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
    st.markdown("Load macromolecules and ligand libraries for pre-processing curation.")
    
    st.divider()
    st.markdown("### 1. Target Macromolecule")
    uploaded_pdb = st.file_uploader("Upload Structure (.pdb)", type=["pdb"], help="Upload any protein or receptor target for structural verification.")
    
    pdb_content = None
    if uploaded_pdb is not None:
        pdb_content = uploaded_pdb.getvalue().decode("utf-8")
    elif os.path.exists("data/1fkn.pdb"):
        with open("data/1fkn.pdb", "r") as f: pdb_content = f.read()
    elif os.path.exists("1fkn.pdb"):
        with open("1fkn.pdb", "r") as f: pdb_content = f.read()
        
    st.divider()
    st.markdown("### 2. Ligand Curation Library")
    input_method = st.radio("Input Method:", ["PubChem Name Search (AI/REST)", "Upload Structure Files (.csv, .sdf, .mol)"])
    
    lib_df = None
    if input_method == "PubChem Name Search (AI/REST)":
        drug_input = st.text_area("Type chemical names (comma separated):", value="Aspirin, Ibuprofen, Donepezil, Caffeine, Imatinib", help="Our engine automatically queries NIH PubChem with typo-correction to fetch 3D chemical structures.")
        if drug_input:
            with st.spinner("Fetching live chemical structures..."):
                lib_df = fetch_from_pubchem(drug_input)
    else:
        uploaded_file = st.file_uploader("Upload Library (.csv, .sdf, .mol, .pdb)", type=["csv", "sdf", "mol", "pdb"])
        if uploaded_file is not None:
            with st.spinner("Parsing chemical file formats..."):
                lib_df = parse_uploaded_ligands(uploaded_file)

    st.divider()
    st.markdown("### 📊 ADME Filter Criteria")
    max_lipinski = st.slider("Max Lipinski Violations", 0, 4, 1)
    filter_druglike = st.checkbox("Require Veber Drug-Likeness", value=True)

# --- MAIN DASHBOARD AREA ---
st.title("🧬 CADD Pre-Processing & Library Curation Workbench")
st.markdown("""
**A lightweight, browser-based gateway for computational drug discovery pipelines.**  
Running 3D grid-search docking (AutoDock Vina) and thermodynamic simulations (GROMACS) requires intensive CPU compute that belongs in command-line automation. This interactive GUI serves as the essential pre-processing step: visually validate 3D receptor structures, standardize raw multi-format chemical libraries (.sdf, .mol), and apply ADME-Tox cheminformatics filters before exporting curated datasets to downstream command-line scripts.
""")

tab1, tab2 = st.tabs(["🏛️ 3D Macromolecule Inspector (SBDD)", "💊 Library Curation & ADME-Tox Engine (LBDD)"])

# --- TAB 1: 3D MOLECULAR VIEWER ---
with tab1:
    st.header("Structure-Based Visual Verification")
    if pdb_content:
        col1, col2 = st.columns([3, 1])
        with col2:
            st.markdown("#### Rendering Controls")
            render_style = st.radio("Receptor Style:", ["cartoon", "surface"], index=0)
            color_scheme = st.selectbox("Coloring Scheme:", ["spectrum", "chain", "secondary structure"], index=0)
            st.markdown("---")
            st.caption("🔬 **Navigation:**\n* **Rotate:** Left-click + drag\n* **Zoom:** Scroll wheel\n* **Active Site:** Green carbon sticks highlight bound heteroatoms/ligands.")
        with col1:
            components.html(render_3d_pdb(pdb_content, style=render_style, color_by=color_scheme.split()[0]), height=540, width=800)
    else:
        st.info("👋 Upload a `.pdb` receptor coordinate file in the sidebar to inspect its 3D architecture, secondary structure, and co-crystallized ligand binding pockets.")

# --- TAB 2: ADME-TOX SCREENING ---
with tab2:
    st.header("Ligand-Based Profiling & Curation")
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
                st.caption("Filtered by Lipinski Rule of 5 and Veber rotational/TPSA thresholds.")
            with col2:
                st.dataframe(filtered_df, use_container_width=True)
                
            csv_export = filtered_df.to_csv(index=False).encode('utf-8')
            st.download_button("📥 Export Curated Screening Report (.csv)", data=csv_export, file_name="curated_ligand_library.csv", mime="text/csv", help="Download this clean CSV to feed directly into command-line docking pipelines.")
    else:
        st.info("👈 Enter drug names in the sidebar or upload a structural library (.sdf, .mol, .csv) to standardize formats and calculate ADME descriptors.")