'''
python v_analysis.py --help
'''

import argparse
import os
import subprocess
import pandas as pd
from Bio.Seq import Seq

def create_fasta(data, output_dir, input_file):
    base_name = os.path.basename(input_file).replace('.tsv', '_igblast.fasta')
    fasta_filename = os.path.join(output_dir, base_name)
    with open(fasta_filename, 'w') as fasta_file:
        for index, entry in data.iterrows():
            header = f"{entry['subject acc.ver']}|{entry['s. start']}|{entry['s. end']}"
            sequence = f"{entry['5prime1k']}{entry['subject_seq']}{entry['3prime50']}"
            fasta_file.write(f">{header}\n{sequence}\n")
    return fasta_filename

def call_igblast(infile, species, output_dir, input_file):
    """ Calls BLAST-based NCBI tool IgBLAST for mapping and annotating gene elements from BCR mRNA sequences."""

    # Define the output file names
    outnametsv = infile.replace('.fasta', '.tsv')
    
    # Define the IgBLAST directory and executable path
    igblast_dir = '/stor/work/Georgiou/Sharing_Folder/IgBLAST_1.21.0/ncbi-igblast-1.21.0'
    igblastn = os.path.join(igblast_dir, 'bin/igblastn')
    script_dir = os.getcwd()
    
    # Construct the command to call IgBLAST with the appropriate parameters and redirect output
    command = (
        f"cd {igblast_dir}; "
        f"{igblastn} "
        f"-germline_db_V internal_data/{species}/{species}_V "
        f"-germline_db_J internal_data/{species}/{species}_J "
        f"-germline_db_D internal_data/{species}/{species}_D "
        f"-auxiliary_data optional_file/{species}_gl.aux "
        f"-organism {species} "
        f"-query {script_dir}/{infile} "
        f"-outfmt 19 "
        f"-out {script_dir}/{outnametsv} "
    )
    
    print('\tCalling IgBLAST to assign FRs & CDRs...')
    subprocess.run(command, shell=True, check=True)
    return outnametsv

def merge_igblast_results(original_data, igblast_output, output_dir):
    # Read the IgBLAST output
    igblast_data = pd.read_csv(igblast_output, sep='\t')
    
    # Ensure the 'sequence_id' used in merging is available in both dataframes
    original_data['sequence_id'] = original_data.apply(
        lambda row: f"{row['subject acc.ver']}|{row['s. start']}|{row['s. end']}", axis=1)
    
    # Merge the dataframes on 'sequence_id'
    relevant_columns = ['sequence_id', 'fwr1', 'fwr1_aa', 'cdr1', 'cdr1_aa', 'fwr2', 'fwr2_aa', 'cdr2', 'cdr2_aa', 'fwr3', 'fwr3_aa', 'cdr3', 'cdr3_aa']
    merged_data = pd.merge(original_data, igblast_data[relevant_columns], on='sequence_id', how='left')
    
    return merged_data

def read_leader_sequences(filepath):
    # Read the TSV file containing the leader sequences
    leader_data = pd.read_csv(filepath, sep='\t', header=0, index_col=None)
    # Convert all sequences to uppercase
    for col in leader_data.columns:
        leader_data[col] = leader_data[col].str.upper()
    return leader_data

def process_leaders(data, leader_sequences, num_best):
    print('\tAligning leader consensus sequences...')
    # Iterate over each row in the DataFrame
    for index, row in data.iterrows():
        gene_type = row['gene_type']
        # Extend 5prime1k with the first 5 nucleotides of subject_seq
        five_prime_sequence = row['5prime1k'] + row['subject_seq'][:5]
        
        # Retrieve input leader consensus sequences, which will be queries
        leader_l1 = leader_sequences[f'{gene_type.lower()}_l1'].iloc[0]
        leader_l2 = leader_sequences[f'{gene_type.lower()}_l2'].iloc[0]

        # Find upstream location & sequences of best alignment to the leader consensus
        best_matches_l1, min_dists_l1 = find_lowest_hamming_window(five_prime_sequence, leader_l1, num_best, False)
        best_matches_l2, min_dists_l2 = find_lowest_hamming_window(five_prime_sequence, leader_l2, num_best, True, len(five_prime_sequence))

        # Update DataFrame with best matches, surrounding sequences, and actual nucleotide matches
        update_dataframe_with_matches(index, data, best_matches_l1, best_matches_l2, min_dists_l1, min_dists_l2)
    return data

def find_lowest_hamming_window(target_sequence, query_sequence, num_best, is_l2=False, sequence_length=0):
    # Crawl along the target sequence and gather all the alignment scores
    distances = [(calculate_hamming_distance(target_sequence[i:i + len(query_sequence)], query_sequence), i)
                 for i in range(len(target_sequence) - len(query_sequence) + 1)]
    
    # Order the alignments by score
    sorted_distances = sorted(distances, key=lambda x: x[0])
    unique_distances = sorted(set([dist for dist, _ in sorted_distances]))[:num_best]

    best_matches = []
    min_distances = []
    for dist in unique_distances:
        match_count = sum(1 for distance, _ in sorted_distances if distance == dist)
        for distance, start in sorted_distances:
            if distance == dist:
                end = start + len(query_sequence)
                if not is_l2 or (is_l2 and (sequence_length - end) <= 55):
                    surrounding_sequences = {
                        'window': target_sequence[start:end],
                        'splice_donor': target_sequence[end:end+10] if end + 10 <= len(target_sequence) else target_sequence[end:],
                        'splice_acceptor': target_sequence[start-10:start] if start >= 10 else target_sequence[:start]
                    }
                    best_matches.append((start, end, surrounding_sequences))
                    if len(best_matches) >= match_count * num_best:
                        break
        # Append each distance match_count times to the list
        min_distances.extend([dist] * match_count)

    return best_matches, min_distances

def calculate_hamming_distance(seq1, seq2):
    """ Calculate the Hamming distance between two sequences of equal length, ignoring positions with 'N'. """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")
    return sum(ch1 != ch2 for ch1, ch2 in zip(seq1, seq2) if 'N' not in (ch1, ch2))

def update_dataframe_with_matches(index, data, best_matches_l1, best_matches_l2, min_dists_l1, min_dists_l2, min_length=38):
    # Extract windows and positions
    L1_windows_positions = [(seqs['window'], start, end) for start, end, seqs in best_matches_l1]
    L2_windows_positions = [(seqs['window'], start, end) for start, end, seqs in best_matches_l2]
    
    # Create concatenated leader nucleotide sequences only where L1 ends before L2 starts
    concatenated_sequences = []
    trimmed_L1s = []
    for l1, _, end_l1 in L1_windows_positions:
        for l2, start_l2, _ in L2_windows_positions:
            if end_l1 < start_l2:  # Ensure L1 ends before L2 starts
                # This adjustment is in case L1 alignment is too long (due to consensus being too long)
                adjusted_L1s = [l1[:i] for i in range(len(l1), min_length - 1, -1)]
                for new_l1 in adjusted_L1s:
                    trimmed_L1s.append(new_l1)
                    concatenated_sequences.append(new_l1 + l2)
    
    # Keep track of the trimmed L1 that gave rise to a given L1+L2 combination
    data.at[index, 'L1_nt_adjusted'] = '; '.join(trimmed_L1s)
    
    # Store concatenated sequences in the dataframe
    data.at[index, 'full_leaders_nt'] = '; '.join(concatenated_sequences)
    
    # Translate full leader sequences
    translated_sequences = [str(Seq(seq).translate()) for seq in concatenated_sequences]
    data.at[index, 'full_leaders_aa'] = '; '.join(translated_sequences)
    
    # Create columns for best L1 & L2 matches (and surrounding sequences)
    data.at[index, 'L1_best_matches'] = '; '.join([f"{start}-{end}" for _, start, end in L1_windows_positions])
    data.at[index, 'L2_best_matches'] = '; '.join([f"{start}-{end}" for start, end, _ in L2_windows_positions])
    data.at[index, 'L1_nt'] = '; '.join([seq for seq, _, _ in L1_windows_positions])
    data.at[index, 'L2_nt'] = '; '.join([seq for seq, _, _ in L2_windows_positions])
    data.at[index, 'L1_splice_donor'] = '; '.join([seqs['splice_donor'] for _, _, seqs in best_matches_l1])
    data.at[index, 'L2_splice_acceptor'] = '; '.join([seqs['splice_acceptor'] for _, _, seqs in best_matches_l2])
    data.at[index, 'L1_start'] = '; '.join([seqs['splice_acceptor'] for _, _, seqs in best_matches_l1])
    data.at[index, 'L1_min_hamming'] = '; '.join(map(str, min_dists_l1))
    data.at[index, 'L2_min_hamming'] = '; '.join(map(str, min_dists_l2))
    
def create_fasta_for_signalp(data, output_dir):
    fasta_path = os.path.join(output_dir, 'leader_candidates.fasta')
    with open(fasta_path, 'w') as fasta_file:
        skipped_count = 0  # Keep track of how many sequences are skipped
        for index, row in data.iterrows():
            sequences = [seq for seq in row['full_leaders_aa'].split('; ') if seq.strip()]  # Ensure no empty or whitespace-only strings
            descriptions = [f">seq{index}_{i}" for i, seq in enumerate(sequences) if seq.strip()]  # Sync descriptions with sequences
            
            for desc, seq in zip(descriptions, sequences):
                if seq.strip():  # Double-check to catch any potential empty entries
                    fasta_file.write(f"{desc}\n{seq}\n")
                else:
                    skipped_count += 1
                    print(f"Skipped empty sequence at index {index}, descriptor {desc}")
                    
        print(f"Total skipped sequences: {skipped_count}")  # Log how many sequences were skipped
    return fasta_path

def run_signalp(fasta_path, output_dir):
    signalp_output_prefix = os.path.join(output_dir, 'leader_candidates')
    cmd = f"signalp -fasta {fasta_path} -org euk -format short -prefix {signalp_output_prefix}"
    subprocess.run(cmd, shell=True, check=True)
    return f"{signalp_output_prefix}_summary.signalp5"

def parse_signalp_output(signalp_output_path):
    signalp_scores = pd.read_csv(signalp_output_path, sep='\t', skiprows=1)
    return signalp_scores[['# ID', 'SP(Sec/SPI)']]

def map_signalp_scores_to_data(data, signalp_scores):
    # Extracting 'index' and 'sub_index' from the '# ID' and assume the ID format in SignalP output is 'seq{index}_{sub_index}'
    signalp_scores['index'] = signalp_scores['# ID'].apply(lambda x: int(x.split('_')[0].strip('seq')))
    signalp_scores['sub_index'] = signalp_scores['# ID'].apply(lambda x: int(x.split('_')[1]))

    # Creating a dictionary with (index, sub_index) as keys and scores as values
    score_dict = signalp_scores.set_index(['index', 'sub_index'])['SP(Sec/SPI)'].to_dict()

    # Applying scores back to data DataFrame by generating scores for each peptide in 'full_leaders_aa'
    data['full_leaders_score'] = data.apply(lambda row: '; '.join(
        [str(score_dict.get((row.name, i), 'NaN')) for i in range(len(row['full_leaders_aa'].split('; ')))])
        , axis=1)

    return data

def find_best_leader_combination(data):
    # Assuming that the scores are already sorted or you compute the max here
    data['best_leader_info'] = data.apply(lambda row: max(zip(
        row['full_leaders_score'].split('; '), 
        row['full_leaders_nt'].split('; '),
        row['full_leaders_aa'].split('; '), 
        row['L1_nt_adjusted'].split('; ')), key=lambda x: float(x[0])), axis=1)
    return data

def main(args):
    output_dir = os.path.dirname(args.input_file)
    # Read input TSV
    data = pd.read_csv(args.input_file, sep='\t')
    # Read leader sequences
    leader_sequences = read_leader_sequences(args.leader_file)
    
    # Perform IgBLAST analysis
    fasta_file = create_fasta(data, output_dir, args.input_file)
    igblast_output = call_igblast(fasta_file, args.species, output_dir, args.input_file)
    # Merge FR & CDR data with the original putative Ig gene loci data
    data = merge_igblast_results(data, igblast_output, output_dir)
    
    # Find candidate leader sequences
    data = process_leaders(data, leader_sequences, args.num_best)
    
    # First round of SignalP
    print('\tCalling SignalP-5.0 to analyze candidate leader sequences combinations..')
    fasta_path = create_fasta_for_signalp(data, output_dir)
    signalp_output = run_signalp(fasta_path, output_dir)
    signalp_scores = parse_signalp_output(signalp_output)
    data = map_signalp_scores_to_data(data, signalp_scores)
    
    # Find the best leader combination
    data = find_best_leader_combination(data)

    # Remove temporary columns before saving the final output (for Excel viewing ease)
    columns_to_drop = ['L1_nt_adjusted', 'full_leaders_nt', 'full_leaders_aa']
    data.drop(columns=columns_to_drop, inplace=True)

    # Save the final annotated data
    final_output_file_path = os.path.join(output_dir, os.path.basename(args.input_file).replace('.tsv', '_annotations.tsv'))
    data.to_csv(final_output_file_path, sep='\t', index=False)
    print(f"Final annotated output available at: {final_output_file_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TSV file of partially annotated BLAST hits for FR, CDR, and leader peptide identification.")
    parser.add_argument('input_file', type=str, help='Input TSV file containing the required columns.')
    parser.add_argument('--species', type=str, default='human', help='Species for IgBLAST reference dataset.')
    parser.add_argument('--leader_file', type=str, help='TSV file containing leader sequences.')
    parser.add_argument('--num_best', type=int, default=3, help='Number of best Hamming distances to consider when aligning putative leader/signal sequences.')
    args = parser.parse_args()
    main(args)

