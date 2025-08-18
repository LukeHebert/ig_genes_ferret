'''
python j_analysis.py J_rss_annotations.tsv
'''

import argparse
import pandas as pd
from Bio.Seq import Seq
import re

def select_index_pair(index_str, column_type):
    """ Select the best index pair based on the type of column ('5p' or '3p'). """
    pairs = index_str.split(';')
    valid_pairs = [pair for pair in pairs if len(pair.split(',')) == 2]  # Validate pairs

    if len(valid_pairs) == 0:
        return None

    if column_type == '5p':
        # Closest to 50 from the end for '5p'
        return min(valid_pairs, key=lambda x: abs(int(x.strip('()').split(', ')[1]) - 50))
    else:
        # Closest to 0 from the start for '3p'
        return min(valid_pairs, key=lambda x: abs(int(x.strip('()').split(', ')[0])))

def create_gene_segment(row):
    five_prime = row['5prime50']
    subject_seq = row['subject_seq']
    three_prime = row['3prime50']

    five_p_idx = select_index_pair(row['combined_5p_indexes'], '5p')
    three_p_idx = select_index_pair(row['combined_3p_indexes'], '3p')

    if not five_p_idx or not three_p_idx:
        return ""

    five_p = int(five_p_idx.strip('()').split(', ')[1])
    gene_segment = "" # gene_segment will be the Jregion + downstream 50 nt (to catch the splice site)
    segment_start = max(0, five_p - len(five_prime))
    if five_p > len(five_prime):
        gene_segment = subject_seq[segment_start:] + three_prime
    else:
        gene_segment = five_prime[five_p:] + subject_seq + three_prime
    return gene_segment

def find_motif_translations(row):
    gene_segment = row['gene_segment']
    gene_type = row['gene_type']
    translations = [str(Seq(gene_segment[i:]).translate(to_stop=False)) for i in range(3)]
    if gene_type in ["IGKJ", "IGLJ"]:
        motif_regex = re.compile("F[GC].G")
    else:
        motif_regex = re.compile("WG.G")
    motif_info = []
    for index, translation in enumerate(translations):
        for match in motif_regex.finditer(translation):
            start_pos = match.start()
            nucleotide_index = index + start_pos * 3
            motif_info.append(f"{index+1}_{start_pos}_{nucleotide_index}")
    return translations, motif_info

def check_splice_site(row, motif_info):
    gene_type = row['gene_type']
    gene_segment = row['gene_segment']
    results = []
    for info in motif_info:
        frame, amino_index, nucleotide_index = map(int, info.split('_'))
        # Adjust offset based on gene_type
        if gene_type in ["IGLJ", "IGKJ"]:
            check_index = nucleotide_index + 19 + 12  # `+ 12` is for the J-MOTIF codons' length
        else:
            check_index = nucleotide_index + 22 + 12  # `+ 12` is for the J-MOTIF codons' length
        if len(gene_segment) > check_index + 2:
            splice_region = gene_segment[max(0, check_index - 2):min(len(gene_segment), check_index + 4)]
            is_gt = gene_segment[check_index:check_index + 2] == "GT"
            results.append((is_gt, splice_region))
        else:
            results.append((False, ""))
    return results

def process_data(file_path):
    df = pd.read_csv(file_path, delimiter='\t')
    df['combined_5p_indexes'] = df['5p_idx_exact'].astype(str) + ';' + df['5p_idx_fuzzy'].astype(str)
    df['combined_3p_indexes'] = df['3p_idx_exact'].astype(str) + ';' + df['3p_idx_fuzzy'].astype(str)

    df['gene_segment'] = df.apply(create_gene_segment, axis=1)
    df['translations'], df['motif_frame_aa_nt'] = zip(*df.apply(find_motif_translations, axis=1))
    results = df.apply(lambda row: check_splice_site(row, row['motif_frame_aa_nt']), axis=1)
    df['splice_d_22'] = results.apply(lambda x: any(item[0] for item in x))
    df['splice_d_22_seq'] = results.apply(lambda x: ' '.join([item[1] for item in x if item[1]]))

    for i in range(1, 4):
        df[f'GSframe{i}'] = df['translations'].apply(lambda x: x[i-1])

    df.drop(columns=['translations', 'combined_5p_indexes', 'combined_3p_indexes', 'rss_score'], inplace=True)

    output_path = file_path.replace('.tsv', '_ORFs.tsv')
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Processed file saved as: {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Process TSV file to find specific motifs and splice sites.")
    parser.add_argument("file_path", type=str, help="Path to the input TSV file.")
    args = parser.parse_args()

    process_data(args.file_path)

if __name__ == "__main__":
    main()

