import os
import shutil

def organize_repository():
    """Restructures a flat workspace into a modular bioinformatics repository."""
    print("[*] Initializing repository cleanup and folder structuring...")
    
    # Define target architecture
    directories = {
        "data": [".pdb", ".fasta", ".sdf", ".mol"],
        "scripts": [".py", ".pml"],
        "results/alignments": [".aln", ".meg"],
        "results/docking": [".pdbqt", "vina_config.txt", "receptor_clean.pdb"],
        "results/md_simulation": [".mdp", ".top", ".gro", ".xvg"],
        "results/reports": [".csv", ".png"]
    }
    
    # Create target directories
    for folder in directories.keys():
        os.makedirs(folder, exist_ok=True)
        
    # Get all files in root directory
    root_files = [f for f in os.listdir(".") if os.path.isfile(f)]
    
    # Files to keep explicitly in the root directory
    root_exceptions = {"organize_workspace.py", "README.md", "environment.yml", ".gitignore"}
    
    moved_count = 0
    for filename in root_files:
        if filename in root_exceptions or filename.startswith("."):
            continue
            
        file_ext = os.path.splitext(filename)[1].lower()
        destination = None
        
        # Check explicit filenames first, then fall back to extensions
        for folder, rules in directories.items():
            if filename in rules or file_ext in rules:
                # Special handling: keep our source protein in data, but clean receptor in docking
                if filename == "1fkn.pdb":
                    destination = "data"
                elif filename == "receptor_clean.pdb":
                    destination = "results/docking"
                elif file_ext == ".pdb":
                    destination = "data"
                else:
                    destination = folder
                break
                
        if destination:
            target_path = os.path.join(destination, filename)
            # Handle potential overwrites gracefully
            if os.path.exists(target_path):
                os.remove(target_path)
            shutil.move(filename, target_path)
            print(f"  [->] Moved: {filename:30} -> {destination}/")
            moved_count += 1
            
    print(f"\n[+] Cleanup complete! Successfully organized {moved_count} artifacts into clean subdirectories.")
    print("[*] Don't forget to move your scripts into 'scripts/' once you're done editing!")

if __name__ == "__main__":
    organize_repository()