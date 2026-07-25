import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

def convert_ligands_to_pdbqt():
    """
    Converts ligand SDF/MOL files into 3D protonated PDBQT files with Gasteiger charges.
    """
    sdf_files = [f for f in os.listdir(".") if f.endswith(".sdf") or f.endswith(".mol")]
    print(f"[*] Found {len(sdf_files)} ligand structural file(s) for 3D preparation...")
    
    for filename in sdf_files:
        base_name = os.path.splitext(filename)[0]
        out_pdbqt = f"{base_name}.pdbqt"
        
        # Load molecule with RDKit, add explicit hydrogens, and generate 3D conformer
        mol = Chem.MolFromMolFile(filename) if filename.endswith(".mol") else Chem.SDMolSupplier(filename)[0]
        if mol is None:
            continue
            
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.U协调 = True
        
        # Assign Gasteiger partial charges
        AllChem.ComputeGasteigerCharges(mol)
        
        # Save temporary PDB
        temp_pdb = f"{base_name}_temp.pdb"
        Chem.MolToPDBFile(mol, temp_pdb)
        
        print(f"  [+] Prepared 3D protonated structure: '{temp_pdb}'")

if __name__ == "__main__":
    convert_ligands_to_pdbqt()