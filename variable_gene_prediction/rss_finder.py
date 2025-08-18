'''
This script takes a single TSV file containing BLASTn hits that have been 
annotated with surrounding sequence data in an attempt to find immunoglobulin
gene features such as RSSs, reading frames, splice sites & 
signal peptide sequences (for V genes). The result is a new set of columns with
these RSS annotations.

python rss_finder.py annotated_hits.tsv rss_seqs.pkl output_name.tsv --num_workers 100
'''
 
 
import argparse
import pandas as pd
import pickle
from Bio.Seq import Seq
import logging
from datetime import datetime
from itertools import product
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed


def setup_logging(output_file):
    # Create a log filename based on the output filename and the current datetime
    log_filename = f"{output_file.replace('.tsv','')}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    # Setup logging to file and console
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    file_handler = logging.FileHandler(log_filename)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(levelname)s - %(message)s'))
    
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(logging.Formatter('%(message)s'))
    
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    
    return log_filename


def load_data(tsv_file, pickle_file):
    df = pd.read_csv(tsv_file, sep='\t')
    df = clean_sequences(df)
    
    with open(pickle_file, 'rb') as f:
        rss_dict = pickle.load(f)
        
    logging.info("BLAST/seq data and RSS patterns loaded successfully.")
    return df, rss_dict


def clean_sequences(df):
    # Fill NaN with empty strings
    cols_to_clean = ['5prime50', '3prime50', 'subject_seq']
    for col in cols_to_clean:
        df[col] = df[col].fillna('')
        df[col] = df[col].astype(str)  # Ensure all data in these columns are strings
    return df
    

def process_sequences(df, rss_dict, num_workers):
    '''Analyzes df rows in parallel by calling process_row() in tandem. Returns
    essentially the original df with new RSS columns added.'''
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Convert DataFrame to a list of dictionaries
        row_dicts = df.to_dict('records')
        
        # Submit jobs with indices
        futures = {executor.submit(process_row, row_dict, rss_dict): i for i, row_dict in enumerate(row_dicts)}
        
        # Collect results in order of submission
        rss_results = [None] * len(row_dicts)  # Pre-initialize a list of the correct size
        for future in as_completed(futures):
            index = futures[future]  # Retrieve the index for this future
            rss_results[index] = future.result()  # Place the result in the corresponding index

    # Convert list of dictionaries to DataFrame
    rss_cols = pd.DataFrame(rss_results)
    
    # Combine the original DataFrame with the new RSS columns
    combined_df = pd.concat([df, rss_cols], axis=1)
    
    # Add a helpful column showing spacer sizes for identified RSSs
    combined_df['5p_exact_spacers'] = combined_df['5p_idx_exact'].apply(which_spacer)
    combined_df['5p_fuzzy_spacers'] = combined_df['5p_idx_fuzzy'].apply(which_spacer)
    combined_df['3p_exact_spacers'] = combined_df['3p_idx_exact'].apply(which_spacer)
    combined_df['3p_fuzzy_spacers'] = combined_df['3p_idx_fuzzy'].apply(which_spacer)

    return combined_df


def process_row(row_dict, rss_dict):
    '''Analyzes a row's sequence columns for presence of recombination 
    signal sequences, by referencing rss_dict. Creates a small dictionary of any
    successfully parsed RSSs, which will become dataframe columns for appending 
    to the original input BLAST hit file's dataframe.'''
    found = {'5p_rss_exact': '', '5p_idx_exact': '', '3p_rss_exact': '', '3p_idx_exact': '',
            '5p_rss_fuzzy': '', '5p_idx_fuzzy': '', '3p_rss_fuzzy': '', '3p_idx_fuzzy': ''}
    subject_seq = row_dict.get('subject_seq', '')
    gene_type = row_dict.get('gene_type', '')

    for neighbor_field in ['5prime50', '3prime50']:
        neighbor_seq = row_dict.get(neighbor_field, None)
        if neighbor_seq is not None:
            direction = neighbor_field[:2]
            combined_seq, offset = combine_sequences(neighbor_seq, subject_seq, direction)

            # First, attempt to find exact matches for all spacer lengths
            exact_found = False
            for spacer_length in ['12', '23']:
                rss, indices = find_exact_matches(combined_seq, rss_dict, spacer_length, offset, direction)
                if rss:
                    if found[f"{direction}_rss_exact"]:
                        found[f"{direction}_rss_exact"] += ';'+str(rss)
                        found[f"{direction}_idx_exact"] += ';'+str(indices)
                    else:
                        found[f"{direction}_rss_exact"] = str(rss)
                        found[f"{direction}_idx_exact"] = str(indices)
                    exact_found = True
            
            # If no exact matches were found, proceed to find fuzzy matches
            if not exact_found:
                for spacer_length in ['12', '23']:
                    rss, indices = find_fuzzy_matches(combined_seq, rss_dict, spacer_length, offset, direction, gene_type)
                    if rss:
                        if found[f"{direction}_rss_fuzzy"]:
                            found[f"{direction}_rss_fuzzy"] += ';'+str(rss)
                            found[f"{direction}_idx_fuzzy"] += ';'+str(indices)
                        else:
                            found[f"{direction}_rss_fuzzy"] = str(rss)
                            found[f"{direction}_idx_fuzzy"] = str(indices)

    return found


def combine_sequences(neighbor_seq, subject_seq, side):
    '''Takes either the 5' or 3' sequence neighboring the alignment subject_seq
    and combines part of the subject_seq to the neighbor_seq in case an RSS 
    slightly overlaps the BLAST alignment sequence'''
    
    # Convert NaN to empty string
    neighbor_seq = '' if pd.isna(neighbor_seq) else neighbor_seq
    subject_seq = '' if pd.isna(subject_seq) else subject_seq
    
    if side == '5p': # (Variable called `direction` elsewhere)
        # Get up to 25bp from the 5' side of the subject sequence
        aligned_piece = subject_seq[:min(25, len(subject_seq))]
        combined_seq = neighbor_seq + aligned_piece  # Add to the start of neighbor_seq
        offset = 0  # No offset needed as aligned_piece is at the end
    elif side == '3p':
        # Get up to 25bp from the 3' side of the subject sequence
        aligned_piece = subject_seq[-min(25, len(subject_seq)):]
        combined_seq = aligned_piece + neighbor_seq  # Add to the end of neighbor_seq
        offset = -len(aligned_piece) # For index adjustment of any located RSS
    return combined_seq, offset


def find_exact_matches(sequence, rss_dict, spacer_length, offset, direction):
    '''Searches for RSS patterns based on the configuration 7 + 12|23 + 9 or its
    reverse. Adjusts search based on the direction (5' or 3'), which determines 
    the order of searching heptamers and nonamers.'''
    sequence = sequence.upper()
    known_heptamers = {seq.upper() for seq in rss_dict[spacer_length]['7']}
    known_nonamers = {seq.upper() for seq in rss_dict[spacer_length]['9']}
    
    # Determine heptamer-nonamer search order based on direction
    if direction == '5p':
        first_set, second_set = known_nonamers, known_heptamers
        first_length, second_length = 9, 7
    else:
        first_set, second_set = known_heptamers, known_nonamers
        first_length, second_length = 7, 9

    # Execute search based on the determined order
    for first_pattern in first_set: # first_patter is a nonamer if looking at 5' flank, heptamer for 3' flank
        first_pattern_rc = str(Seq(first_pattern).reverse_complement())
        # The tuple with a boolean helps keep track of pattern strand so that e.g. a revcomp nonamer does not get paired with a forward heptamer
        for search_first, search_first_rc in [(first_pattern, False), (first_pattern_rc, True)]:
            pos_first = sequence.find(search_first)
            if pos_first != -1:
                second_start = pos_first + len(search_first) + int(spacer_length)
                second_end = second_start + second_length
                if second_end <= len(sequence):
                    maybe_second = sequence[second_start:second_end]
                    for second in second_set:
                        second_rc = str(Seq(second).reverse_complement())
                        for second_pattern, is_rc in [(second, False), (second_rc, True)]:
                            if maybe_second == second_pattern and search_first_rc == is_rc:
                                return (f"{search_first}-{sequence[pos_first+len(search_first):second_start]}-{maybe_second}",
                                        (pos_first + offset, second_end + offset))
    return None, None


def find_fuzzy_matches(sequence, rss_dict, spacer_length, offset, direction, gene_type):
    '''Searches for RSS patterns based on the configuration 7 + 12|23 + 9 or its
    reverse. Adjusts search based on the direction (5' or 3'), which determines 
    the order of searching heptamers and nonamers. Allows for flexible/inexact 
    matches to the reference RSS sequences in order to anticipate undocumented 
    inter-species variations.'''
    sequence = sequence.upper()
    allowed_mismatches = {'7': 1, '9': 1} if gene_type == 'IGHD' else {'7': 2, '9': 3}

    known_heptamers = {seq.upper() for seq in rss_dict[spacer_length]['7']}
    known_nonamers = {seq.upper() for seq in rss_dict[spacer_length]['9']}

    # Determine heptamer-nonamer search order based on direction
    if direction == '5p':
        first_set, second_set = known_nonamers, known_heptamers
        first_length, second_length = 9, 7
    else:
        first_set, second_set = known_heptamers, known_nonamers
        first_length, second_length = 7, 9

    # Fuzzy search based on the determined order
    for first_pattern in first_set: # first_patter is a nonamer if looking at 5' flank, heptamer for 3' flank
        first_pattern_rc = str(Seq(first_pattern).reverse_complement())
        # The tuple with a boolean helps keep track of pattern strand so that e.g. a revcomp nonamer does not get paired with a forward heptamer
        for search_first, search_first_rc in [(first_pattern, False), (first_pattern_rc, True)]:
            for start_pos in range(len(sequence) - int(first_length) + 1):
                test_first = sequence[start_pos:start_pos + int(first_length)]
                if hamming_distance(test_first, search_first) <= allowed_mismatches[str(first_length)]:
                    second_start = start_pos + int(first_length) + int(spacer_length)
                    second_end = second_start + int(second_length)
                    if second_end <= len(sequence):
                        test_second = sequence[second_start:second_end]
                        for second in second_set:
                            second_rc = str(Seq(second).reverse_complement())
                            for second_pattern, is_rc in [(second, False), (second_rc, True)]:
                                if hamming_distance(test_second, second_pattern) <= allowed_mismatches[str(second_length)] and search_first_rc == is_rc:
                                    return (f"{search_first}-{sequence[start_pos+len(search_first):second_start]}-{test_second}",
                                            (start_pos + offset, second_end + offset))
    return None, None


def hamming_distance(s1, s2):
    '''Calculate the Hamming distance between two strings of equal length'''
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def which_spacer(entry):
    '''Assumes `entry` is a semicolon-separated string of pairs like 
    '(-14, 14);(-14, 25)' & returns the absolute values of the differences of
    index pairs. `entry` can have an arbitrary number of index pairs.'''
    if not entry:
        return ''
    
    pairs = entry.replace(" ", "").split(";")
    differences = []
    
    for pair in pairs:
        if pair.startswith('(') and pair.endswith(')'):
            numbers = pair[1:-1].split(",")
            if len(numbers) == 2 and all(num.isdigit() or (num.startswith('-') and num[1:].isdigit()) for num in numbers):
                diff = abs(int(numbers[0]) - int(numbers[1]))
                differences.append(str(diff))
            else:
                differences.append('')  # Append an empty string or some placeholder for invalid data
        else:
            differences.append('')  # Append an empty string or some placeholder for invalid data
    
    return ";".join(differences)


def score_rss(row):
    '''This function helps score RSSs based on 1) whether they were detected and
    2) how far away they are from a BLAST alignment (functional RSSs are directly
    adjacent)'''
    # Check the primary category
    categories = {
        1: ['5p_rss_exact', '5p_idx_exact', '3p_rss_exact', '3p_idx_exact', '5p_rss_fuzzy', '5p_idx_fuzzy', '3p_rss_fuzzy', '3p_idx_fuzzy'],
        2: ['5p_rss_exact', '5p_idx_exact', '3p_rss_exact', '3p_idx_exact'],
        3: ['5p_rss_fuzzy', '5p_idx_fuzzy', '3p_rss_fuzzy', '3p_idx_fuzzy'],
        4: ['5p_rss_exact', '5p_idx_exact', '3p_rss_exact', '3p_idx_exact', '5p_rss_fuzzy', '5p_idx_fuzzy', '3p_rss_fuzzy', '3p_idx_fuzzy'],
    }
    
    primary_score = 5  # default to the lowest category
    for score, cols in categories.items():
        if all(pd.notna(row[col]) and row[col] != '' for col in cols):
            primary_score = score
            break
    
    # Calculate secondary score for exact columns
    min_diff_5p = min((abs(int(b) - 50) for a, b in (pair.strip('()').split(',') for pair in str(row['5p_idx_exact']).split(';') if pair)), default=0)
    min_diff_3p = min((abs(int(a) - 0) for a, b in (pair.strip('()').split(',') for pair in str(row['3p_idx_exact']).split(';') if pair)), default=0)
    
    secondary_score = min_diff_5p + min_diff_3p
    
    # Combine scores: primary score dominates, secondary sorts within primary
    return primary_score * 1000 - secondary_score


def check_file_suffix(filename):
    if not filename.endswith('.tsv'):
        raise argparse.ArgumentTypeError("Output file must have a '.tsv' suffix.")
    return filename


def main(args):
    # Set up log file
    log_filename = setup_logging(args.output_file)
    logging.info(f"Command-line arguments: {sys.argv}")
    
    # Load input datafiles
    df, rss_dict = load_data(args.tsv_file, args.pickle_file)
    
    # Identify any RSS-like patterns in sequences that flank BLAST hits
    processed_df = process_sequences(df, rss_dict, args.num_workers)

    # Sort the results by gene type and by "realness" of the identified RSSs
    gene_type_order = ['IGHV', 'IGHD', 'IGHJ', 'IGKV', 'IGKJ', 'IGLV', 'IGLJ']
    gene_type_categorical = pd.Categorical(processed_df['gene_type'], categories=gene_type_order, ordered=True)
    processed_df['gene_type'] = gene_type_categorical
    # Calculate the distance scores
    processed_df['rss_score'] = processed_df.apply(score_rss, axis=1)
    # Calculate minimum of 's. start' and 's. end' for locus sorting
    processed_df['min_start_end'] = processed_df[['s. start', 's. end']].min(axis=1)
    # Finally sort the output
    processed_df.sort_values(by=['gene_type', 'rss_score', 'evalue'], ascending=[True,False,True], inplace=True)

    # Log the count of rows per gene_type before filtering
    pre_filter_gene_type_counts = processed_df['gene_type'].value_counts()
    for gene_type, count in pre_filter_gene_type_counts.items():
        logging.info(f"{gene_type}: {count} rows before filtering")

    # Save outputs separated by gene type
    for gene_type in gene_type_order:
        gene_specific_df = processed_df[processed_df['gene_type'] == gene_type]
        gene_specific_output_file = f"{args.output_file.replace('.tsv', '')}_{gene_type}.tsv"
        gene_specific_df.to_csv(gene_specific_output_file, sep='\t', index=False)
        logging.info(f"Unfiltered output for {gene_type} saved to {gene_specific_output_file}. Contains all original rows of this gene_type, annotated.")

    # Filter the DataFrame to include rows with at least one RSS found,
    # unless the row's gene_type is IGHD, in which case two RSSs must be found
    
    # Individual conditions
    has_5p_exact = processed_df['5p_rss_exact'].str.strip() != ''
    has_5p_fuzzy = processed_df['5p_rss_fuzzy'].str.strip() != ''
    has_3p_exact = processed_df['3p_rss_exact'].str.strip() != ''
    has_3p_fuzzy = processed_df['3p_rss_fuzzy'].str.strip() != ''

    # Grouped conditions (for readability)
    has_any_5p = has_5p_exact | has_5p_fuzzy
    has_any_3p = has_3p_exact | has_3p_fuzzy
    has_both_5p_3p = has_any_5p & has_any_3p

    # Condition for non-IGHD genes
    condition_non_ighd = (processed_df['gene_type'] != 'IGHD') & (has_any_5p | has_any_3p)

    # Condition for IGHD genes
    condition_ighd = (processed_df['gene_type'] == 'IGHD') & has_both_5p_3p

    # Filter the DataFrame
    rss_filtered_df = processed_df[condition_non_ighd | condition_ighd]
    
    # Log the count of rows per gene_type in the filtered dataframe
    gene_type_counts = rss_filtered_df['gene_type'].value_counts()
    for gene_type, count in gene_type_counts.items():
        logging.info(f"{gene_type}: {count} rows found with requisite RSS(s)")
    # Save filtered outputs separated by gene type
    for gene_type in gene_type_order:
        filtered_gene_specific_df = rss_filtered_df[rss_filtered_df['gene_type'] == gene_type]
        filtered_output_file = f"{args.output_file.replace('.tsv', '')}_{gene_type}_FoundOnly.tsv"
        filtered_gene_specific_df.to_csv(filtered_output_file, sep='\t', index=False)
        logging.info(f"Filtered output for {gene_type} saved to {filtered_output_file}. Contains only rows with the RSS(s) necessary for this gene type.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Search for RSS patterns adjacent to immunoglobulin loci.")
    parser.add_argument("tsv_file", help="Path to the input TSV file containing BLASTn results and adjacent sequences.")
    parser.add_argument("pickle_file", help="Path to the binary pickle file containing RSS patterns of the format {'12':{'7':set(),'9':set()},'23':{'7':set(),'9':set()}}.")
    parser.add_argument("output_file", help="Path to the output TSV file containing additional columns for RSS patterns. Must have a '.tsv' suffix.", type=check_file_suffix)
    parser.add_argument("--num_workers", type=int, default=4, help="Number of worker processes to use for parallel processing.")
    args = parser.parse_args()
    main(args)

