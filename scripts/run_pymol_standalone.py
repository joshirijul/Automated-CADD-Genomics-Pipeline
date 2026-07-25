import os
import sys
import shutil
import subprocess

def locate_pymol_executable():
    """Scans Windows environment PATH and common installation directories for PyMol executable."""
    # 1. Check if it's already registered in Windows PATH under common alias names
    for alias in ["pymol", "PyMOLWin", "pymol.exe", "PyMOLWin.exe"]:
        found = shutil.which(alias)
        if found:
            return found

    # 2. Scan standard Windows installation directories
    user_home = os.path.expanduser("~")
    common_paths = [
        r"C:\Program Files\PyMOL\PyMOLWin.exe",
        r"C:\Program Files\PyMOL\pymol.exe",
        r"C:\Program Files (x86)\PyMOL\PyMOLWin.exe",
        r"C:\Program Files (x86)\PyMOL\pymol.exe",
        os.path.join(user_home, r"AppData\Local\Programs\PyMOL\PyMOLWin.exe"),
        os.path.join(user_home, r"AppData\Local\Programs\PyMOL\pymol.exe"),
        os.path.join(user_home, r"anaconda3\envs\cadd_env\Library\bin\pymol.bat"),
        os.path.join(user_home, r"miniconda3\envs\cadd_env\Library\bin\pymol.bat"),
    ]

    for path in common_paths:
        if os.path.exists(path):
            return path
            
    return None

def execute_standalone_mapping(pml_script="map_active_site.pml", output_png="bace1_active_site_topography.png"):
    if not os.path.exists("1fkn.pdb"):
        print("[!] Error: '1fkn.pdb' not found in your workspace.")
        return

    print("[*] Locating standalone PyMol binary on Windows...")
    pymol_exe = locate_pymol_executable()
    
    if not pymol_exe:
        print("\n[!] Could not automatically locate PyMol binary.")
        print("[*] FALLBACK MANUAL STEP:")
        print("    1. Open your PyMol application normally.")
        print(f"    2. In the PyMol upper command line, type: cd {os.path.abspath('.')}")
        print(f"    3. Then type: @{pml_script}")
        return

    print(f"[+] Found PyMol executable: '{pymol_exe}'")
    print(f"[*] Dispatching headless background rendering process for '{pml_script}'...")
    
    try:
        # Launch PyMol as a separate background process, passing the -c (headless) flag and macro file
        cmd_args = [pymol_exe, "-c", pml_script]
        result = subprocess.run(cmd_args, capture_output=True, text=True, check=True)
        
        if os.path.exists(output_png):
            print(f"\n[+] Success! Active site topography rendered to '{output_png}'")
        else:
            print("[!] Subprocess finished, but PNG was not detected. Console output:")
            print(result.stdout)
            
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Subprocess failed with return code {e.returncode}:")
        print(e.stderr or e.stdout)

if __name__ == "__main__":
    execute_standalone_mapping()