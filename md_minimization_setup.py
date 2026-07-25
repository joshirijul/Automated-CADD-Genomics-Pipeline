import os
import matplotlib.pyplot as plt
import seaborn as sns

def generate_mdp_files():
    """Generates standard GROMACS parameter files for ionization and minimization."""
    
    ions_mdp = """t; ions.mdp - used as input into grompp to generate ions.tpr
integrator  = steep         ; Algorithm (steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when Fmax < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size (nm)
nsteps      = 50000         ; Maximum number of steps
nstlist         = 1         ; Frequency to update neighbor list
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
"""

    minim_mdp = """; minim.mdp - parameters for steepest-descents energy minimization
integrator  = steep         ; Algorithm (steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when Fmax < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size (nm)
nsteps      = 50000         ; Maximum number of minimization steps to perform
nstlist         = 1         ; Frequency to update neighbor list
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list
coulombtype     = PME       ; Particle Mesh Ewald for long-range electrostatics
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions
"""

    with open("ions.mdp", "w") as f:
        f.write(ions_mdp)
    with open("minim.mdp", "w") as f:
        f.write(minim_mdp)
        
    print("[+] Successfully generated 'ions.mdp' and 'minim.mdp' in workspace.")

def plot_minimization_curve(xvg_file="potential.xvg", output_png="energy_minimization_curve.png"):
    """Parses GROMACS energy XVG output to plot potential energy relaxation."""
    if not os.path.exists(xvg_file):
        print(f"[*] '{xvg_file}' not detected yet.")
        print("    -> Execute the GROMACS terminal pipeline first, then re-run this script to plot results!")
        return

    print(f"[*] Parsing GROMACS energy file '{xvg_file}'...")
    steps = []
    energies = []
    
    with open(xvg_file, 'r') as f:
        for line in f:
            # Ignore XVG headers and metadata lines
            if line.startswith("#") or line.startswith("@"):
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                steps.append(float(parts[0]))
                energies.append(float(parts[1]))
                
    if not steps:
        print("[!] Error: Could not extract numerical data from XVG file.")
        return

    print(f"[+] Extracted {len(steps)} energy frames.")
    
    plt.figure(figsize=(8, 6))
    sns.set_style("whitegrid")
    plt.plot(steps, energies, color="#1f77b4", linewidth=2, label="Potential Energy ($E_{pot}$)")
    
    plt.title("GROMACS Steepest-Descents Energy Minimization — BACE1", fontsize=14, fontweight="bold")
    plt.xlabel("Minimization Step", fontsize=12)
    plt.ylabel("Potential Energy ($E_{pot}$ in kJ/mol)", fontsize=12)
    plt.legend(loc="upper right", frameon=True)
    plt.tight_layout()
    plt.savefig(output_png, dpi=300)
    print(f"[+] Publication-ready minimization curve saved to '{output_png}'\n")

if __name__ == "__main__":
    generate_mdp_files()
    plot_minimization_curve()