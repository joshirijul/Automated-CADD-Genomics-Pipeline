import os
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PPBuilder

def generate_ramachandran(pdb_file="1fkn.pdb", output_png="ramachandran_1fkn.png"):
    if not os.path.exists(pdb_file):
        print(f"[!] Could not find '{pdb_file}'. Ensure your target PDB is in the workspace.")
        return

    print(f"[*] Parsing structure and calculating dihedral angles for '{pdb_file}'...")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("BACE1", pdb_file)
    
    phi_angles = []
    psi_angles = []
    
    # Use BioPython's Polypeptide Builder to extract backbone geometry
    ppb = PPBuilder()
    for model in structure:
        for chain in model:
            for pp in ppb.build_peptides(chain):
                phi_psi = pp.get_phi_psi_list()
                for phi, psi in phi_psi:
                    if phi is not None and psi is not None:
                        # Convert radians to degrees
                        phi_angles.append(np.degrees(phi))
                        psi_angles.append(np.degrees(psi))
                        
    print(f"[+] Successfully extracted {len(phi_angles)} residue conformations.")
    
    # Plotting the Ramachandran space
    plt.figure(figsize=(8, 8))
    plt.scatter(phi_angles, psi_angles, s=15, c="#1f77b4", alpha=0.7, edgecolors="none", label="Residues")
    
    # Define reference lines and quadrants
    plt.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    plt.axvline(0, color="gray", linestyle="--", linewidth=0.8)
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    
    plt.title("Ramachandran Plot — Human BACE1 (1FKN)", fontsize=14, fontweight="bold")
    plt.xlabel("Phi ($\phi$) Angle (degrees)", fontsize=12)
    plt.ylabel("Psi ($\psi$) Angle (degrees)", fontsize=12)
    plt.grid(True, linestyle=":", alpha=0.6)
    
    # Add a visual bounding box for core beta-sheet and alpha-helix regions
    plt.axvspan(-180, -30, ymin=0.5, ymax=1.0, color="green", alpha=0.08, label="Beta-Sheet Region")
    plt.axvspan(-180, -30, ymin=0.0, ymax=0.5, color="red", alpha=0.08, label="Alpha-Helix Region")
    
    plt.legend(loc="upper right")
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"[+] Ramachandran plot saved as publication-ready graphic: '{output_png}'\n")

if __name__ == "__main__":
    generate_ramachandran()