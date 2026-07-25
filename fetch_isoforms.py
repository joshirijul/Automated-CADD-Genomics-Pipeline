import urllib.request
import os

def fetch_bace1_isoforms(output_file='bace1_isoforms.fasta'):
    print(f"[*] Querying Uniprot for human BACE1 (P56817) isoforms...")

    url= "https://rest.uniprot.org/uniprotkb/search?query=accession:P56817&format=fasta&includeIsoform=true"

    try:
        urllib.request.urlretrieve(url, output_file)

        with open(output_file, 'r') as f:
            content = f.read()
            sequence_count = content.count('>')

        print(f"[+] Success ! Downloaded {sequence_count} isoform sequences.")
        print(f"[+] Sequences saved to '{os.path.abspath(output_file)}'")
    except Exception as e:
        print(f"[ERROR] Sequence retrieval failed: {str(e)}")

if __name__ == "__main__":
    fetch_bace1_isoforms()
        
