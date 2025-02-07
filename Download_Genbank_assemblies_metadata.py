import argparse
import os
import csv
from Bio import Entrez, SeqIO
from tqdm import tqdm

def read_accessions(file_path):
    with open(file_path, "r") as f:
        return [line.strip() for line in f.readlines()[1:] if line.strip()]

def fetch_metadata(accession_numbers):
    metadata = []
    for acc in tqdm(accession_numbers, desc="Fetching metadata", unit="accession"):
        try:
            handle = Entrez.efetch(db="nucleotide", id=acc, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            meta_info = {
                "Accession": record.id,
                "Organism": record.annotations.get("organism", "N/A"),
                "Definition": record.description,
                "Length": len(record.seq),
                "Source": record.annotations.get("source", "N/A"),
                "Collection Date": record.features[0].qualifiers.get("collection_date", ["N/A"])[0] if record.features else "N/A",
                "Isolation Source": record.features[0].qualifiers.get("isolation_source", ["N/A"])[0] if record.features else "N/A",
                "Geo Location": record.features[0].qualifiers.get("geo_loc_name", ["N/A"])[0] if record.features else "N/A"
            }
            metadata.append(meta_info)
        except Exception as e:
            print(f"Error fetching {acc}: {e}")
    return metadata

def save_metadata(metadata, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "metadata.csv")
    
    with open(output_file, "w", newline="") as csvfile:
        fieldnames = ["Accession", "Organism", "Definition", "Length", "Source", "Collection Date", "Isolation Source", "Geo Location"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(metadata)
    
    print(f"Metadata saved to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Fetch metadata for a list of GenBank accession numbers.")
    parser.add_argument("-t", "--accessions", required=True, help="Text file containing GenBank accession numbers with a header 'Accessions'.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for the metadata CSV file.")
    
    args = parser.parse_args()
    Entrez.email = "Benjamin.guinet95@gmail.com"  # Required by NCBI Entrez API
    
    accession_numbers = read_accessions(args.accessions)
    metadata = fetch_metadata(accession_numbers)
    save_metadata(metadata, args.output)

if __name__ == "__main__":
    main()
