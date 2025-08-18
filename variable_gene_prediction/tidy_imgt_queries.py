'''
Takes a directory as input. The directory should contain an arbitrary number of 
FASTA files all named with a common prefix e.g. "IGHV_dog.fasta" & 
"IGHV_human.fasta". This script concatenates the input FASTAs into a single file
and removes any `.` charactrers in the input FASTA file record sequences.
'''


import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def clean_sequence(sequence):
    """Remove all period characters from a sequence."""
    return sequence.replace('.', '')

def combine_fasta(directory, output_file):
    """Combine all fasta files in the directory into one, cleaning sequences along the way."""
    # Initialize a list to hold all cleaned SeqRecord objects
    all_records = []

    # Process each file
    for fasta_file in os.listdir(directory):
        if fasta_file.endswith('.fasta') or fasta_file.endswith('.fa'):
            file_path = os.path.join(directory, fasta_file)
            # Read sequences from each file
            for record in SeqIO.parse(file_path, "fasta"):
                # Clean the sequence
                cleaned_sequence = clean_sequence(str(record.seq))
                # Create a new SeqRecord with the cleaned sequence
                new_record = SeqRecord(Seq(cleaned_sequence), id=record.id, description=record.description)
                all_records.append(new_record)


    # Write all cleaned SeqRecords to the output file
    SeqIO.write(all_records, output_file, "fasta")
    print(f"Combined and cleaned fasta written to {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Combine fasta files into a single file with cleaned sequences.")
    parser.add_argument("directory", type=str, help="Directory containing fasta files.")
    
    args = parser.parse_args()

    # Assume all fasta files share the same four-letter prefix
    fasta_files = [f for f in os.listdir(args.directory) if f.endswith('.fasta')]
    if not fasta_files:
        print("No fasta files found in the directory.")
        return
    
    prefix = fasta_files[0][:4]
    output_file = os.path.join(args.directory, f"{prefix}_all.fasta")

    combine_fasta(args.directory, output_file)

if __name__ == "__main__":
    main()
