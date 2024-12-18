import argparse
import pandas as pd

def consolidate_columns(df):
    imgt_cols = [col for col in df.columns if col.startswith('imgt_')]
    ref_cols = [col for col in df.columns if col.startswith('ref_')]
    
    # Map each ref column to its corresponding imgt column
    for ref_col in ref_cols:
        imgt_col = 'imgt_' + ref_col[4:]  # Corresponding imgt column name
        if imgt_col in imgt_cols:
            # Consolidate columns if both ref and imgt versions exist
            df[imgt_col] = df[imgt_col].fillna(df[ref_col])
        else:
            # Rename ref columns to imgt if imgt version doesn't exist
            df.rename(columns={ref_col: imgt_col}, inplace=True)
    
    # Now remove all original ref columns after processing
    ref_cols_to_remove = [col for col in df.columns if col.startswith('ref_')]
    df.drop(columns=ref_cols_to_remove, inplace=True, errors='ignore')

    return df


def remove_unwanted_columns(df):
    # Remove unwanted columns by name and pattern matching
    columns_to_drop = ['matched', 'key', 'matches', 'seq_gapped'] + [col for col in df.columns if 'tata_box' in col or 'octamer' in col or col.startswith('blast_')]
    df.drop(columns=columns_to_drop, inplace=True, errors='ignore')
    return df

def preserve_columns(df):
    if 'query acc.ver' in df.columns:
        col_index = df.columns.get_loc('query acc.ver') + 1
        always_preserve = df.columns[:col_index].tolist()
    else:
        always_preserve = []

    additional_columns = {
        'V': ['best_leader_info', '5prime1k', '3prime50', '3p_rss_exact', '3p_idx_exact', '3p_rss_fuzzy', '3p_idx_fuzzy', '3p_fuzzy_spacers', 'fwr1','fwr1_aa','cdr1','cdr1_aa','fwr2','fwr2_aa','cdr2','cdr2_aa','fwr3','fwr3_aa','cdr3','cdr3_aa'],
        'D': ['frame1', 'frame2', 'frame3', '5prime50', '5p_rss_exact', '5p_idx_exact', '5p_exact_spacers', '5p_rss_fuzzy', '5p_idx_fuzzy', '3prime50', '3p_rss_exact', '3p_idx_exact', '3p_rss_fuzzy', '3p_idx_fuzzy', '3p_fuzzy_spacers'],
        'J': ['5prime50', '5p_rss_exact', '5p_idx_exact', '5p_exact_spacers', '5p_rss_fuzzy', '5p_idx_fuzzy', 'motif_frame_aa_nt', 'splice_donor_distance', 'splice_donor_seq', 'gene_segment']
    }

    # Strip any spaces that might exist in the gene_type column
    df['gene_type'] = df['gene_type'].str.strip()
    # A direct row by row check approach
    for key, cols in additional_columns.items():
        mask = df['gene_type'].str.contains(key, na=False)
        always_preserve += [col for col in cols if col in df.columns and mask.any()]

    return df.filter(items=always_preserve)


def rename_and_reformat(df):
    # Rename columns and reformat values as specified
    rename_dict = {
        'GSframe1': 'frame1',
        'GSframe2': 'frame2',
        'GSframe3': 'frame3',
        'splice_d_22': 'splice_donor_distance',
        'splice_d_22_seq': 'splice_donor_seq'
    }
    df.rename(columns=rename_dict, inplace=True)
    df['splice_donor_distance'] = df.apply(lambda row: 22 if row['gene_type'] == 'IGHJ' else 19 if row['gene_type'] in ['IGKJ', 'IGLJ'] else 'NA', axis=1)
    return df


def process_files(file_list):
    processed_dataframes = []
    for file in file_list:
        df = pd.read_csv(file, sep='\t')
        df = consolidate_columns(df)
        df = remove_unwanted_columns(df)
        df = rename_and_reformat(df)
        df = preserve_columns(df)
        df.reset_index(drop=True, inplace=True)
        processed_dataframes.append(df)
    
    concatenated_df = pd.concat(processed_dataframes, ignore_index=True)
    return concatenated_df

def main():
    parser = argparse.ArgumentParser(description="Process TSV files to format and consolidate data.")
    parser.add_argument("files", nargs='+', help="List of TSV files to process")
    args = parser.parse_args()
    final_df = process_files(args.files)
    final_df.to_csv('all_overlapping.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
