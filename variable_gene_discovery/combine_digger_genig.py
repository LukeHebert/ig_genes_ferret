'''
Takes:
    1 digger annotated gene loci .csv output file
    1 "genig" (Hebert et al Seeger) annotated gene loci .tsv output file
Does:
    finds how many loci overlap & how many were unique to each input file based on scaffold & start/end positions
    for those loci that overlap, merges & saves that data as a new .tsv file
    for those loci that are unique to genig/Hebert-Seeger output, create a subset of that .tsv file as a new .tsv file
Output:
    .txt file with overlap/exclusion counts
    .tsv with merged overlapping data
'''

import argparse
import pandas as pd

def clean_subject_acc_ver(value):
    """ Remove periods and convert to lowercase. """
    return value.replace('.', '').lower()

def check_overlap(start1, end1, start2, end2):
    """ Check if two ranges overlap, allowing for ranges where start could be greater than end. """
    start1, end1 = sorted([start1, end1])
    start2, end2 = sorted([start2, end2])
    return max(start1, start2) <= min(end1, end2)

def main(digger_file, genig_file, gene_type):
    # Load data
    digger_df = pd.read_csv(digger_file)
    genig_df = pd.read_csv(genig_file, sep='\t')

    # Convert gene_type to uppercase and filter data
    gene_type = gene_type.upper()
    digger_df = digger_df[digger_df['gene_type'] == gene_type]
    genig_df = genig_df[genig_df['gene_type'] == gene_type]

    # Clean the GenIg 'subject acc.ver' column
    genig_df['subject acc.ver'] = genig_df['subject acc.ver'].apply(clean_subject_acc_ver)

    # Prepare for merging and comparison
    digger_df['key'] = digger_df['contig'].str.lower()
    genig_df['key'] = genig_df['subject acc.ver']

    # Initialize a column to track matched rows
    digger_df['matched'] = False

    # Find overlapping entries
    overlap_data = []
    genig_only_data = genig_df.copy()

    for index, dig_row in digger_df.iterrows():
        matched = False
        for _, gen_row in genig_df.iterrows():
            if dig_row['key'] == gen_row['key'] and check_overlap(dig_row['start'], dig_row['end'], gen_row['s. start'], gen_row['s. end']):
                overlap_data.append(pd.concat([dig_row, gen_row], axis=0))
                matched = True
                genig_only_data = genig_only_data.drop(gen_row.name)
        digger_df.at[index, 'matched'] = matched

    # Save the overlapping data
    if overlap_data:
        overlap_df = pd.concat(overlap_data, axis=1).transpose()
        overlap_file = digger_file.replace('.csv', f'_{gene_type}overlap.tsv')
        overlap_df.to_csv(overlap_file, sep='\t', index=False)
        print(f"Overlapping data saved to {overlap_file}")

    # Save GenIg only data
    if not genig_only_data.empty:
        genig_only_file = genig_file.replace('.tsv', '_onlyGenIg.tsv')
        genig_only_data.to_csv(genig_only_file, sep='\t', index=False)
        print(f"GenIg only data saved to {genig_only_file}")

    # Save Digger only data
    digger_only_df = digger_df[digger_df['matched'] == False]
    if not digger_only_df.empty:
        digger_only_file = digger_file.replace('.csv', '_onlyDigger.tsv')
        digger_only_df.to_csv(digger_only_file, sep='\t', index=False)
        print(f"Digger only data saved to {digger_only_file}")

    # Save summary
    summary = {
        'Digger not in GenIg': len(digger_only_df),
        'GenIg not in Digger': len(genig_only_data),
        'Overlapping': len(overlap_data)
    }

    summary_file = f"{gene_type}_overlap_summary.txt"
    with open(summary_file, 'w') as f:
        for key, value in summary.items():
            f.write(f"{key}: {value}\n")
    print(f"Summary saved to {summary_file}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process and find overlaps between Digger and GenIg datasets.')
    parser.add_argument('digger_file', type=str, help='The CSV file from the Digger tool.')
    parser.add_argument('genig_file', type=str, help='The TSV file from the GenIg tool.')
    parser.add_argument('gene_type', type=str, help='Gene type to filter the data on (will be converted to uppercase).')
    args = parser.parse_args()

    main(args.digger_file, args.genig_file, args.gene_type)


