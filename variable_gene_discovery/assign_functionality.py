'''
Guided by figure 7 of 2010 IMGT LIG-Motif publication (doi: 10.1186/1471-2105-11-223): 
V
    no stop codons for L-PART1 & V-EXON (check if either Digger or Genig found L1 & L2)
    are splice sites found, & do they combine L-PART1 + V-EXON in-frame?
    Are 1st-CYS, conserved-TRP, 2nd-CYS in the same reading frame? (requires IMGT-VQUEST)
    Is the RSS spacer length correct for its gene_type/chain?
    Are the heptamer & nonamer identical to previously reported ones in functional genes? (requires RSS pickle file)
D
    at least 1 open reading frame? (requires GenIg reading frames columns)
    is the RSS spacer lengths correct?
    have the RSSs been previously reported around functional genes? (requires RSS pickle file)
J
    splice site found?
    no stop codons in the motif-containing reading frame?
    splice site in frame?
    is the conserved AA J-motif present?
    is the RSS spacer length correct?
    have the RSSs been previously reported around functional genes?
'''

import argparse
import pandas as pd
import numpy as np
import pickle

def parse_leader_info(info):
    if info:
        return info.strip("()").strip("'").split(', ')
    return []

def load_nested_dict(file_path):
    with open(file_path, 'rb') as file:
        return pickle.load(file)
   
def assign_v_functionality(row, known_rss_seqs):
    row['v_no_stops'] = False
    row['v_splices'] = False
    row['v_conserved'] = False
    row['rss_spacer'] = False
    row['rss_known'] = False
    row['final_functionality'] = "Pseudogene"  # Default value
       
    # Process coding sequences
    if row['aa']:
        row['coding_aa'] = row['aa']
    else:
        # Gathering sequences from various columns
        leader_info = parse_leader_info(row['best_leader_info']) if pd.notna(row['best_leader_info']) else []
        coding_seq = leader_info + [row[c] for c in ['fwr1_aa', 'cdr1_aa', 'fwr2_aa', 'cdr2_aa', 'fwr3_aa', 'cdr3_aa'] if pd.notna(row[c])]
        if coding_seq:
            row['coding_aa'] = ''.join(coding_seq)
        else:
            row['coding_aa'] = np.nan  # Ensure NaN is explicitly set if no sequences are available

    # v_no_stops check - Ensure coding_aa is not NaN before checking for '*'
    if pd.notna(row['coding_aa']) and '*' not in row['coding_aa']:
        row['v_no_stops'] = True
    else:
        row['v_no_stops'] = False
    
    # v_splices check
    if pd.notna(row['l_part1']) and pd.notna(row['l_part2']):
        row['v_splices'] = True

    # v_conserved check with gene_type specific criteria
    if row['functional'] == 'Functional':
        row['v_conserved'] = True
    else:
        # Check for absence of certain phrases in 'notes'
        notes_issues = any(phrase in row['notes'] for phrase in ['First cysteine not found', 'Conserved Trp not found', 'Second cysteine not found'])
        if not notes_issues:
            row['v_conserved'] = True
        else:
            # Additional checks on amino acid positions depending on gene_type
            conserved = True  # Assume conserved unless proven otherwise
            
            if 'IGHV' in row['gene_type']:
                if pd.notna(row['fwr1_aa']) and len(row['fwr1_aa']) >= 22:
                    conserved &= row['fwr1_aa'][21] == 'C'  # 1-based index 22 for IGHV
                if pd.notna(row['fwr2_aa']) and len(row['fwr2_aa']) >= 3:
                    conserved &= row['fwr2_aa'][2] == 'W'  # 1-based index 3 for all
                if pd.notna(row['fwr3_aa']) and len(row['fwr3_aa']) >= 38:
                    conserved &= row['fwr3_aa'][37] == 'C'  # 1-based index 38 for IGHV
            
            elif 'IGKV' in row['gene_type']:
                if pd.notna(row['fwr1_aa']) and len(row['fwr1_aa']) >= 24:
                    conserved &= row['fwr1_aa'][22] == 'C' or row['fwr1_aa'][23] == 'C'  # 1-based index 23 or 24 for IGKV
                if pd.notna(row['fwr2_aa']) and len(row['fwr2_aa']) >= 3:
                    conserved &= row['fwr2_aa'][2] == 'W'  # Same for all
                if pd.notna(row['fwr3_aa']) and len(row['fwr3_aa']) >= 36:
                    conserved &= row['fwr3_aa'][35] == 'C' or row['fwr3_aa'][32] == 'C'  # 1-based index 36 or 33 for IGKV
            
            elif 'IGLV' in row['gene_type']:
                if pd.notna(row['fwr1_aa']) and len(row['fwr1_aa']) >= 22:
                    conserved &= row['fwr1_aa'][21] == 'C'  # 1-based index 22 for IGLV
                if pd.notna(row['fwr2_aa']) and len(row['fwr2_aa']) >= 3:
                    conserved &= row['fwr2_aa'][2] == 'W'  # Same for all
                if pd.notna(row['fwr3_aa']) and len(row['fwr3_aa']) >= 38:
                    conserved &= row['fwr3_aa'][35] == 'C' or row['fwr3_aa'][34] == 'C' or row['fwr3_aa'][37] == 'C'  # 1-based index 36, 35, or 38 for IGLV

            row['v_conserved'] = conserved

    # rss_spacer calculation
    digger_spacer_len = None
    inhouse_spacer_len = None
    if pd.notna(row['3_rss_start']) and pd.notna(row['3_rss_end']):
        digger_spacer_len = abs(row['3_rss_start'] - row['3_rss_end']) - 15
    elif pd.notna(row['3p_rss_exact']):
        inhouse_spacer_len = len(row['3p_rss_exact'].split('-')[1])
    elif pd.notna(row['3p_rss_fuzzy']):
        inhouse_spacer_len = len(row['3p_rss_fuzzy'].split('-')[1])

    spacer_expected = 23 if row['gene_type'] in ['IGHV', 'IGLV'] else 12 if row['gene_type'] == 'IGKV' else None
    if (digger_spacer_len is not None and digger_spacer_len == spacer_expected) or (inhouse_spacer_len is not None and inhouse_spacer_len == spacer_expected):
        row['rss_spacer'] = True
    else:
        row['rss_spacer'] = False

    # rss_known check
    heptamer_set = known_rss_seqs['heptamer_set']
    nonamer_set = known_rss_seqs['nonamer_set']
    if row['v_heptamer'].upper() in heptamer_set and row['v_nonamer'].upper() in nonamer_set:
        row['rss_known'] = True
    else:
        row['rss_known'] = False

    # Determine final functionality based on conditions
    if all([row[col] for col in ['v_no_stops', 'v_splices', 'v_conserved', 'rss_spacer', 'rss_known']]):
        row['final_functionality'] = "Functional"
    elif all([row[col] for col in ['v_no_stops', 'v_splices']]) and not all([row[col] for col in ['rss_known', 'rss_spacer', 'v_conserved']]):
        row['final_functionality'] = "ORF"
    else:
        row['final_functionality'] = "Pseudogene"

    return row

def assign_d_functionality(row, known_rss_seqs):
    row['d_an_orf'] = False
    row['5_rss_spacer'] = False
    row['3_rss_spacer'] = False
    row['both_rss_spacers'] = False
    row['5_rss_known'] = False
    row['3_rss_known'] = False
    row['both_rss_known'] = False
    
    # d_an_orf
    row['d_an_orf'] = any('*' not in row[col] for col in ['frame1', 'frame2', 'frame3'] if pd.notna(row[col]))

    # 5p rss spacer
    digger_spacer_len = None
    inhouse_spacer_len = None
    if pd.notna(row['5_rss_start']) and pd.notna(row['5_rss_end']):
        digger_spacer_len = abs(row['5_rss_start'] - row['5_rss_end']) - 15
    if pd.notna(row['5p_rss_exact']):
        inhouse_spacer_len = len(row['5p_rss_exact'].split('-')[1])
    elif pd.notna(row['5p_rss_fuzzy']):
        inhouse_spacer_len = len(row['5p_rss_fuzzy'].split('-')[1])

    spacer_expected = 12
    if (digger_spacer_len is not None and digger_spacer_len == spacer_expected) or (inhouse_spacer_len is not None and inhouse_spacer_len == spacer_expected):
        row['5_rss_spacer'] = True
    else:
        row['5_rss_spacer'] = False

    digger_spacer_len = None
    inhouse_spacer_len = None
    # 3p rss spacer
    if pd.notna(row['3_rss_start']) and pd.notna(row['3_rss_end']):
        digger_spacer_len = abs(row['3_rss_start'] - row['3_rss_end']) - 15
    if pd.notna(row['3p_rss_exact']):
        inhouse_spacer_len = len(row['3p_rss_exact'].split('-')[1])
    elif pd.notna(row['3p_rss_fuzzy']):
        inhouse_spacer_len = len(row['3p_rss_fuzzy'].split('-')[1])

    spacer_expected = 12
    if (digger_spacer_len is not None and digger_spacer_len == spacer_expected) or (inhouse_spacer_len is not None and inhouse_spacer_len == spacer_expected):
        row['3_rss_spacer'] = True
    else:
        row['3_rss_spacer'] = False

    # both_rss_spacers
    row['both_rss_spacers'] = (row['5_rss_spacer'] and row['3_rss_spacer'])

    # 3_rss_known
    heptamer_set = known_rss_seqs['heptamer_set']
    nonamer_set = known_rss_seqs['nonamer_set']
    if row['d_5_heptamer'].upper() in heptamer_set and row['d_5_nonamer'].upper() in nonamer_set:
        row['5_rss_known'] = True
    else:
        row['5_rss_known'] = False

    # 3_rss_known
    if row['d_3_heptamer'].upper() in heptamer_set and row['d_3_nonamer'].upper() in nonamer_set:
        row['3_rss_known'] = True
    else:
        row['3_rss_known'] = False

    # both_rss_known
    row['both_rss_known'] = (row['5_rss_known'] and row['3_rss_known'])

    # Determine final functionality based on conditions
    if all([row[col] for col in ['d_an_orf', 'both_rss_spacers', 'both_rss_known']]):
        row['final_functionality'] = "Functional"
    elif row['d_an_orf'] and not all([row[col] for col in ['both_rss_spacers', 'both_rss_known']]):
        row['final_functionality'] = "ORF"
    else:
        row['final_functionality'] = "Pseudogene"

    return row

def assign_j_functionality(row, known_rss_seqs):
    # j_splice
    row['j_splice'] = (row['splice_donor_seq'][2:4] == 'GT') if pd.notna(row['splice_donor_seq']) else False

    # j_no_stops
    row['j_no_stops'] = ('*' not in row['aa']) if pd.notna(row['aa']) else False

    # j_splice_frame
    if row['gene_type'] in ['IGHJ']:
        row['j_splice_frame'] = (row['splice_donor_distance'] == 22)
    else:  # Applies to IGKV and IGLV
        row['j_splice_frame'] = (row['splice_donor_distance'] == 19) if pd.notna(row['splice_donor_distance']) else False

    # j_motif
    row['j_motif'] = (row['motif_frame_aa_nt'] != '') if pd.notna(row['motif_frame_aa_nt']) else False

    
    # rss_spacer calculation
    digger_spacer_len = None
    inhouse_spacer_len = None
    if pd.notna(row['5_rss_start']) and pd.notna(row['5_rss_end']):
        digger_spacer_len = abs(row['5_rss_start'] - row['5_rss_end']) - 15
    if pd.notna(row['5p_rss_exact']):
        inhouse_spacer_len = len(row['5p_rss_exact'].split('-')[1])
    elif pd.notna(row['5p_rss_fuzzy']):
        inhouse_spacer_len = len(row['5p_rss_fuzzy'].split('-')[1])
    spacer_expected = 23 if row['gene_type'] in ['IGHJ', 'IGKJ'] else 12 if row['gene_type'] == 'IGLJ' else None
    if (digger_spacer_len is not None and digger_spacer_len == spacer_expected) or (inhouse_spacer_len is not None and inhouse_spacer_len == spacer_expected):
        row['rss_spacer'] = True
    else:
        row['rss_spacer'] = False
        
    # rss_known check
    heptamer_set = known_rss_seqs['heptamer_set']
    nonamer_set = known_rss_seqs['nonamer_set']
    if row['j_heptamer'].upper() in heptamer_set and row['j_nonamer'].upper() in nonamer_set:
        row['rss_known'] = True
    else:
        row['rss_known'] = False

    # Determine final functionality based on conditions
    if all([row[col] for col in ['j_splice', 'j_no_stops', 'j_splice_frame', 'j_motif', 'rss_spacer', 'rss_known']]):
        row['final_functionality'] = "Functional"
    elif all([row[col] for col in ['j_splice', 'j_no_stops']]) and not all([row[col] for col in ['rss_known', 'rss_spacer', 'j_motif']]):
        row['final_functionality'] = "ORF"
    else:
        row['final_functionality'] = "Pseudogene"

    return row
    
def process_row(row, known_rss_seqs):
    # Determine the gene type and assign functionality accordingly
    if 'V' in row['gene_type']:
        return assign_v_functionality(row, known_rss_seqs)
    elif 'D' in row['gene_type']:
        return assign_d_functionality(row, known_rss_seqs)
    elif 'J' in row['gene_type']:
        return assign_j_functionality(row, known_rss_seqs)
    return row

def read_rss(tsv_file):
    # Load the TSV file using pandas
    df = pd.read_csv(tsv_file, delimiter='\t')

    # Filter and create sets for heptamers and nonamers based on the 'seq_type' column
    heptamer_set = set(df[df['seq_type'] == 'HEPTAMER']['seq'].str.upper())
    nonamer_set = set(df[df['seq_type'] == 'NONAMER']['seq'].str.upper())

    # Create the dictionary with all heptamers and nonamers
    known_rss_seqs = {
        'heptamer_set': heptamer_set,
        'nonamer_set': nonamer_set
    }

    return known_rss_seqs


def sort_cols_rows(df):
    # Define the desired order of columns as specified
    ordered_column_names = [
        "contig", "sense", "start", "end", "start_rev", "end_rev", "gene_type",
        "imgt_match", "imgt_nt_diffs", "imgt_score", "v-gene_aligned_aa", "functional",
        "notes", "seq", "aa", "gene_start", "gene_end", "gene_start_rev", "gene_end_rev",
        "gene_seq", "l_part1_start", "l_part1_start_rev", "l_part1_end", "l_part1_end_rev",
        "l_part1", "l_part2_start", "l_part2_start_rev", "l_part2_end", "l_part2_end_rev",
        "l_part2", "3_rss_start", "3_rss_start_rev", "3_rss_end", "3_rss_end_rev", 
        "5_rss_start", "5_rss_start_rev", "5_rss_end", "5_rss_end_rev", "v_heptamer",
        "v_nonamer", "d_5_heptamer", "d_5_nonamer", "d_3_heptamer", "d_3_nonamer", "j_heptamer",
        "j_nonamer", "likelihood", "evalue", "query acc.ver", "5prime1k", "best_leader_info",
        "fwr1", "fwr1_aa", "cdr1", "cdr1_aa", "fwr2", "fwr2_aa", "cdr2", "cdr2_aa",
        "fwr3", "fwr3_aa", "cdr3", "cdr3_aa", "3prime50", "3p_idx_exact", "3p_rss_exact",
        "3p_idx_fuzzy", "3p_rss_fuzzy", "5prime50", "5p_idx_exact", "5p_rss_exact",
        "5p_idx_fuzzy", "5p_rss_fuzzy", "frame1", "frame2", "frame3", "motif_frame_aa_nt",
        "final_functionality"
    ]

    # Find any columns in df that are not explicitly mentioned in the ordered list
    additional_columns = [col for col in df.columns if col not in ordered_column_names]

    # Combine the ordered columns with any additional columns found in the DataFrame
    new_column_order = ordered_column_names + additional_columns

    # Reorder the DataFrame according to the new column order
    # Ensure only to use column names that are actually present in the DataFrame to avoid KeyErrors
    new_column_order = [col for col in new_column_order if col in df.columns]

    df = df[new_column_order]
    
    # Also sort the rows
    df = df.sort_values(by=['gene_type', 'final_functionality'], ascending=[True, True])

    return df

def main():
    parser = argparse.ArgumentParser(description="Assign functionality to immunoglobulin gene candidates using Digger and GenIg annotations.")
    parser.add_argument("candidates_file", type=str, help="Path to the input TSV file full of Digger- & GenIg- annotated BLAST hits.")
    parser.add_argument("rss_file", type=str, help="Path to the input RSS data TSV file containing all IMGT-reported recombination signal sequences heptamers & nonamers")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.candidates_file, delimiter='\t')
    known_rss_seqs = read_rss(args.rss_file)

    # Apply functionality assignments based on gene_type
    df = df.apply(lambda row: process_row(row, known_rss_seqs), axis=1)
    
    # Order columns
    df = sort_cols_rows(df)
    
    # Save processed data
    df.to_csv(args.candidates_file.replace('.tsv','_FXN_assigned.tsv'), index=False, sep='\t')

if __name__ == "__main__":
    main()
