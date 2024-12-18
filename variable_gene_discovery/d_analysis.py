'''
python d_analysis.py D_rss_annotations.tsv
'''

import argparse
import pandas as pd
from Bio.Seq import Seq
import os

def select_index_pair(index_str, column_type):
    """ Select the best RSS index pair based on the type of column ('5p' or '3p'). """
    pairs = index_str.split(';')
    if len(pairs) == 0:
        return None

    valid_pairs = [pair for pair in pairs if len(pair.split(',')) == 2]  # Ensure each pair has two elements

    if len(valid_pairs) == 0:
        return None

    if column_type == '5p':
        # Choosing index with second integer closest to 50
        return min(valid_pairs, key=lambda x: abs(int(x.strip('()').split(', ')[1]) - 50))
    else:
        # Choosing index with first integer closest to 0
        return min(valid_pairs, key=lambda x: abs(int(x.strip('()').split(', ')[0])))


def create_gene_segment(row):
    """ Generate the putative "gene segment" sequence from the alignment and flanking sequences. """
    five_prime = row['5prime50']
    subject_seq = row['subject_seq']
    three_prime = row['3prime50']

    # Use the select_index_pair function to choose the best indices
    five_p_idx = select_index_pair(row['combined_5p_indexes'], '5p')
    three_p_idx = select_index_pair(row['combined_3p_indexes'], '3p')

    if five_p_idx is None or three_p_idx is None:
        return ""

    # Split to get integers
    five_p = int(five_p_idx.strip('()').split(', ')[1])
    three_p = int(three_p_idx.strip('()').split(', ')[0])

    # Adjusting the indexes as described
    segment_start = max(0, five_p - len(five_prime))
    segment_end = min(len(three_prime), three_p)
    
    # Construct gene_segment based on adjusted indexes
    gene_segment = ""
    if five_p > len(five_prime):
        gene_segment = subject_seq[segment_start:]
    else:
        gene_segment = five_prime[five_p:] + subject_seq

    if three_p < 0:
        gene_segment = gene_segment[:segment_end]
    else:
        gene_segment += three_prime[:three_p]
    
    return gene_segment

def translate_frames(gene_segment):
    """ Translate the gene segment in three frames. """
    if len(gene_segment) % 3 != 0:
        gene_segment = gene_segment[:-(len(gene_segment) % 3)]  # Trimming to make divisible by 3
    
    seq_obj = Seq(gene_segment)
    return [str(seq_obj[i:].translate(to_stop=False)) for i in range(3)]

def process_data(file_path):
    """ Load the TSV, process data and save a new file. """
    df = pd.read_csv(file_path, delimiter='\t')

    # Check for the '28' condition in the specified columns
    condition = df.apply(lambda x: '28' in (str(x['5p_exact_spacers']) + str(x['5p_fuzzy_spacers']) + str(x['3p_exact_spacers']) + str(x['3p_fuzzy_spacers'])), axis=1)
    df = df[condition]

    # Combine indexes into one column for processing
    df['combined_5p_indexes'] = df['5p_idx_exact'].astype(str) + ';' + df['5p_idx_fuzzy'].astype(str)
    df['combined_3p_indexes'] = df['3p_idx_exact'].astype(str) + ';' + df['3p_idx_fuzzy'].astype(str)

    # Apply functions to create new columns
    df['gene_segment'] = df.apply(create_gene_segment, axis=1)
    df[['GSframe1', 'GSframe2', 'GSframe3']] = df['gene_segment'].apply(lambda x: pd.Series(translate_frames(x)))

    # Save output
    output_path = os.path.splitext(file_path)[0] + '_ORFs.tsv'
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Processed file saved as: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Process TSV file to generate gene segments and their reading frames.")
    parser.add_argument("file_path", type=str, help="Path to the input TSV file.")
    args = parser.parse_args()

    process_data(args.file_path)

if __name__ == "__main__":
    main()



