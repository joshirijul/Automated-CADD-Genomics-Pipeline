import os
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select

class ReceptorSelect(Select):
    """Filters out waters and heteroatoms (ligands, ions) to retain only the protein receptor."""
    def accept_residue(self, residue):
        # Keep standard amino acid residues, discard HOH (water) and HETATM ligands
        return residue.id[0] == ' '

class LigandSelect(Select):
    """Isolates non-water heteroatoms (the co-crystallized BACE1 inhibitor)."""
    def accept_residue(self, residue):
        return residue.id[0] != ' ' and residue.get_resname() != 'HOH'

def prepare_receptor_and_grid(pdb_file="1fkn.pdb", config_file="vina_config.txt"):
    if not os.path.exists(pdb_file):
        print(f"[!] Error: '{pdb_file}' not found in workspace.")
        return

    print(f"[*] Parsing structure '{pdb_file}' using BioPython...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("BACE1", pdb_file)
    
    # 1. Clean Receptor Export
    io = PDBIO()
    io.set_structure(structure)
    receptor_pdb = "receptor_clean.pdb"
    io.save(receptor_pdb, ReceptorSelect())
    print(f"[+] Clean protein receptor written to '{receptor_pdb}'")

    # 2. Extract Bound Inhibitor and Compute Center of Geometry
    ligand_coords = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ' and residue.get_resname() != 'HOH':
                    for atom in residue:
                        ligand_coords.append(atom.get_coord())

    if not ligand_coords:
        print("[!] Warning: No bound ligand found. Defaulting to BACE1 catalytic dyad center.")
        # Default BACE1 catalytic dyad (Asp32/Asp228) approximate center
        center_x, center_y, center_z = 16.5, 34.2, 14.8
    else:
        coords_matrix = np.array(ligand_coords)
        center_x, center_y, center_z = np.mean(coords_matrix, axis=0)
        
    center_x, center_y, center_z = round(center_x, 3), round(center_y, 3), round(center_z, 3)
    
    print("\n--- Active Site Grid Box Parameters ---")
    print(f"  Center Coordinates: X = {center_x}, Y = {center_y}, Z = {center_z}")
    print(f"  Box Dimensions:     size_x = 22.0 Å, size_y = 22.0 Å, size_z = 22.0 Å")
    print(f"  Exhaustiveness:     8 (Standard accuracy/speed balance)")

    # 3. Write AutoDock Vina Configuration File
    config_content = f"""# AutoDock Vina Configuration File — BACE1 (1FKN)
receptor = receptor_clean.pdbqt
ligand = donepezil.pdbqt

out = docking_pose_out.pdbqt
log = docking_log.txt

center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = 22.0
size_y = 22.0
size_z = 22.0

exhaustiveness = 8
num_modes = 9
energy_range = 3
"""
    
    with open(config_file, "w") as f:
        f.write(config_content)
        
    print(f"[+] Vina configuration successfully saved to '{config_file}'\n")
    return center_x, center_y, center_z

if __name__ == "__main__":
    prepare_receptor_and_grid()