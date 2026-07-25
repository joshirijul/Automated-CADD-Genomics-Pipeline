import os

def fix_for_mega(input_aln="bace1_muscle.aln", output_meg="bace1_clean.meg"):
    if not os.path.exists(input_aln):
        print(f"[!] Error: Cannot find '{input_aln}'. Ensure you are in the LAB WORK folder.")
        return

    print(f"[*] Reading raw MUSCLE alignment from '{input_aln}'...")
    with open(input_aln, 'r') as f:
        lines = f.readlines()

    with open(output_meg, 'w') as out:
        # Explicitly force MEGA 12 to recognize this as Protein data
        out.write("#Mega\n")
        out.write("!Title BACE1 Isoform Alignment;\n")
        out.write("!Format DataType=Protein;\n\n")

        for line in lines:
            if line.startswith(">"):
                # Strip spaces and symbols (| and -) that crash MEGA's parser
                header = line.strip()[1:]  # Remove the FASTA '>'
                parts = header.split()
                # Create a clean ID: e.g., 'sp|P56817-5|BACE1_HUMAN' -> 'P56817_5_BACE1_HUMAN'
                clean_id = parts[0].replace("|", "_").replace("-", "_")
                out.write(f"\n#{clean_id}\n")
            elif line.strip() and not line.startswith("#") and not line.startswith("!"):
                # Write clean amino acid sequence lines
                out.write(line.strip() + "\n")

    print(f"[+] Success! Generated strict MEGA protein file: '{output_meg}'")

if __name__ == "__main__":
    fix_for_mega()