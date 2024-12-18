'''
Takes BLAST output and uses it to search within a genome assembly FASTA to 
retrieve flanking sequences and annotate the BLAST results with those sequences.
The BLAST results are parsed to determine each hit's immunoglobulin gene type,
which influences which flanking sequences are retrieved (e.g. 50 bp on the 3'
end of a J gene hit would not be useful because J genes do not have 3' RSSs, so
only a 5' 50 bp flanking sequence is retrieved for hits where the query is 
identified as a J gene).

EXAMPLE
python flanking_seqs.py blast_format.tsv genome.fasta
'''

import pandas as pd
from Bio import SeqIO
import argparse

def segment_type(acc_ver):
    # List of gene type identifiers
    gene_types = ["IGHV", "IGKV", "IGLV", "IGHJ", "IGKJ", "IGLJ", "IGHD"]
    # Check if any identifier is in the accession version
    for gene_type in gene_types:
        if gene_type in acc_ver:
            return gene_type
    return "Unknown"  # Default type if no identifier matches

def extract_sequences(tsv_file, fasta_file):
    # Load the TSV file
    df = pd.read_csv(tsv_file, sep='\t')

    # Create empty columns for gene_type, 5prime50, 5prime1k, subject_seq, and 3prime50
    df['gene_type'] = df['query acc.ver'].apply(segment_type)
    df['5prime50'] = ''
    df['5prime1k'] = ''
    df['subject_seq'] = ''
    df['3prime50'] = ''

    # Load the FASTA file
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    # Iterate through the dataframe
    for index, row in df.iterrows():
        print(f'\rProcessing TSV row {index + 1}/{len(df)}...', end='')

        gene_type = row['gene_type']
        subject_acc_ver = row['subject acc.ver']
        start, end = int(row['s. start'])-1, int(row['s. end'])-1
        strand = row['orientation']

        if subject_acc_ver in fasta_sequences:
            seq_record = fasta_sequences[subject_acc_ver].seq
            
            if strand == 'REV':
                seq_record = seq_record.reverse_complement()
                seq_length = len(seq_record)
                # Correctly recalculating 'start' and 'end'
                new_start = seq_length - end - 1
                new_end = seq_length - start - 1
                
                # Ensure the start is always less than the end
                if new_start > new_end:
                    new_start, new_end = new_end, new_start
                
                # Update the start and end with the corrected values
                start, end = new_start, new_end


            # Extract subject_seq
            subject_seq = seq_record[start:end+1]
            df.at[index, 'subject_seq'] = str(subject_seq)

            # Extract 5prime50 and 3prime50 based on gene_type conditions
            if gene_type.endswith(('D', 'J')):
                five_prime_start = max(0, start - 50)
                five_prime_seq = seq_record[five_prime_start:start]
                df.at[index, '5prime50'] = str(five_prime_seq)

            if gene_type.endswith(('V', 'D', 'J')):
                three_prime_end = min(len(seq_record), end + 50)
                three_prime_seq = seq_record[end+1:three_prime_end+1]
                df.at[index, '3prime50'] = str(three_prime_seq)

            # Extract 5prime1k for gene types ending with 'V'
            if gene_type.endswith('V'):
                five_prime1k_start = max(0, start - 1000)
                five_prime1k_seq = seq_record[five_prime1k_start:start]
                df.at[index, '5prime1k'] = str(five_prime1k_seq)

        else:
            print(f'\r\tSkipping accession: {subject_acc_ver}. No matching record was found.')

    # Clear the carriage return line
    print('\r' + ' ' * 50 + '\r', end='')
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Searches a genome FASTA file for BLAST hit locations and annotates the BLAST data with surrounding sequences relevant to immunoglobulin gene analysis.')
    parser.add_argument('tsv_file', type=str, help='Input BLAST result (TSV) file path')
    parser.add_argument('fasta_file', type=str, help='Input FASTA file path')

    args = parser.parse_args()

    df = extract_sequences(args.tsv_file, args.fasta_file)

    outfile_name = args.tsv_file.replace('_distinct.tsv', '_seqs.tsv')
    df.to_csv(outfile_name, sep='\t', index=False)
