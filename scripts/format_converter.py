import os
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors

def convert_workspace_ligands(input_dir='.', output_csv='ligand_smiles_library.csv'):

    supported_extensions=("*.sdf", "*.mol", "*.pdb")
    file_paths = []

    for ext in supported_extensions:
        file_paths.extend(glob.glob(os.path.join(input_dir, ext)))
        file_paths.extend(glob.glob(os.path.join(input_dir, ext.upper())))

    file_paths=sorted(list(set(file_paths)))
    if not file_paths:
        print(f"[!] No .sdf, .mol, or .pdb files found in '{os.path.abspath(input_dir)}'.")
        return None
        
    records = []
    print(f"[*] Scanning active workspace: '{os.path.abspath(input_dir)}'")
    print(f"[*] Found {len(file_paths)} total structure files...")
    
    for path in file_paths:
        filename = os.path.basename(path)
        ext = filename.split(".")[-1].lower()
        
        # MACROMOLECULE FILTER: Skip files > 200 KB (Proteins are typically 500KB+, ligands < 20KB)
        if os.path.getsize(path) > 200 * 1024:
            print(f"  [!] Skipping {filename} (File size > 200KB — detected as protein target)")
            continue
            
        mols = []
        try:
            if ext == "sdf":
                supplier = Chem.SDMolSupplier(path, sanitize=True)
                mols = [m for m in supplier if m is not None]
            elif ext == "mol":
                mol = Chem.MolFromMolFile(path, sanitize=True)
                if mol: mols.append(mol)
            elif ext == "pdb":
                mol = Chem.MolFromPDBFile(path, sanitize=True)
                if mol: mols.append(mol)
                
            if not mols:
                print(f"  [WARNING] Could not parse {filename}. Check valence/stereochemistry.")
                continue
                
            for idx, mol in enumerate(mols):
                mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"{filename.split('.')[0]}_{idx+1}"
                heavy_atoms = mol.GetNumHeavyAtoms()
                
                # SECONDARY FILTER: Skip peptide chains or large fragments > 100 heavy atoms
                if heavy_atoms > 100:
                    print(f"  [!] Skipping {mol_name} in {filename} ({heavy_atoms} heavy atoms — exceeds small-molecule threshold)")
                    continue

                smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
                mw = round(Descriptors.MolWt(mol), 2)
                
                records.append({
                    "File_Name": filename,
                    "Molecule_Name": mol_name,
                    "Canonical_SMILES": smiles,
                    "Heavy_Atoms": heavy_atoms,
                    "Molecular_Weight": mw
                })
                print(f"  [->] Converted: {filename:15} -> {smiles[:35]}...")
                
        except Exception as e:
            print(f"  [ERROR] Failed on {filename}: {str(e)}")
            
    if not records:
        print("\n[!] No valid small molecules were converted. Ensure your ligand files are in this folder.")
        return None
        
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"\n[+] Success! Exported {len(df)} ligand records to '{output_csv}'.")
    return df

if __name__ == "__main__":
    # Scans whichever folder is currently open in your VS Code terminal
    convert_workspace_ligands(input_dir=".")
        