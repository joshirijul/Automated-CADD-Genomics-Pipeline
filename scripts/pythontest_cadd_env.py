import sys
print(f"Testing environment on Python version: {sys.version.split()[0]}\n" + "-"*50)

# 1. Test RDKit Chemoinformatics Pipeline
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    
    # Test SMILES parsing (Aspirin)
    smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Failed to parse SMILES string into Mol object.")
    
    # Test 2D coordinate generation and descriptor calculation
    AllChem.Compute2DCoords(mol)
    mol_wt = Descriptors.MolWt(mol)
    mol_block = Chem.MolToMolBlock(mol)
    
    print(f"[SUCCESS] RDKit Operational")
    print(f"          -> Parsed Aspirin | MolWt: {mol_wt:.2f} g/mol")
    print(f"          -> Generated MOL block ({len(mol_block.splitlines())} lines)")
except Exception as e:
    print(f"[ERROR] RDKit verification failed: {e}")

print("-" * 50)

# 2. Test Biopython Sequence & Molecular Biology Pipeline
try:
    import Bio
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    
    # Test DNA sequence manipulation and central dogma translation
    dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    mrna_seq = dna_seq.transcribe()
    protein_seq = dna_seq.translate(to_stop=True)
    
    record = SeqRecord(dna_seq, id="dry_lab_test", description="Verification sequence")
    
    print(f"[SUCCESS] Biopython v{Bio.__version__} Operational")
    print(f"          -> DNA Record ID: {record.id}")
    print(f"          -> Translated Functional Peptide: {protein_seq}")
except Exception as e:
    print(f"[ERROR] Biopython verification failed: {e}")

print("-" * 50 + "\nVerification complete. Dry lab environment is ready!")