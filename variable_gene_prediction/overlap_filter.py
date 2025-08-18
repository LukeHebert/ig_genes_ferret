'''
1 Parses the input string representing the scientific notation.
2 Assigns the mappings as FWD or REV based on mapped locus orientation
3 Implements a sweep-line algorithm to group overlapping sequence hits
4 Keeps only one row from each overlap group (the longest subject sequence)
5 Exports the resulting dataframes to TSV files
6 Prints updates and create a log file

Example use:
python overlap_filter.py input.tsv
'''


import pandas as pd
import sys
import os

def filter_dataframe(df, evalue):
    return df[df['evalue'] <= evalue]

def add_orientation_column(df):
    # This function adds an 'orientation' column based on s. start and s. end
    df['orientation'] = df.apply(lambda x: 'FWD' if x['s. start'] < x['s. end'] else 'REV', axis=1)
    return df

def determine_intervals_order(df):
    if not df.empty and df.iloc[0]['s. start'] > df.iloc[0]['s. end']:
        return ['s. end', 's. start']
    else:
        return ['s. start', 's. end']

def by_subject_seq(df):
    """
    Groups the dataframe by 'subject acc.ver' and 'orientation', then applies overlap grouping
    within each group using the sweep_line_group() function.
    """
    final_df = pd.DataFrame()

    for _, group_df in df.groupby(['subject acc.ver', 'orientation']):
        processed_group = sweep_line_group(group_df)
        final_df = pd.concat([final_df, processed_group], ignore_index=True)

    return final_df

def sweep_line_group(df):
    '''
    Groups rows in the dataframe based on overlapping ranges in 's. start' and 's. end' columns,
    and selects the row with the largest range from each group. This function uses a sweep-line
    algorithm to efficiently find overlapping intervals.
    '''
    if df.empty:
        return df

    interval_order = determine_intervals_order(df)
    intervals = df[interval_order].to_numpy()
    events = []
    for i, (start, end) in enumerate(intervals):
        events.append((start, True, i))  # Start event
        events.append((end, False, i))   # End event
    events.sort()

    active = set()
    groups = {}
    for position, is_start, index in events:
        if is_start:
            if not active:  # New group starts if no active intervals
                group_key = (index,)
                groups[group_key] = [index]
            else:  # Update existing groups
                for key in list(groups):
                    if set(key).intersection(active):
                        new_key = tuple(sorted(set(key).union({index})))
                        groups[new_key] = groups.pop(key) + [index]
            active.add(index)
        else:
            active.remove(index)

    # Select the row with the largest range for each group
    selected_rows = []
    for group in groups.values():
        max_range_row = max(group, key=lambda i: intervals[i][1] - intervals[i][0])
        selected_rows.append(max_range_row)

    return df.iloc[selected_rows]

def create_log_file(log_path, messages):
    with open(log_path, 'w') as log_file:
        for message in messages:
            log_file.write(message + '\n')

def main():
    if len(sys.argv) != 2:
        print("Usage: python overlap_filter.py <tsv_file>")
        sys.exit(1)

    tsv_file = sys.argv[1]

    print("Reading BLAST file...")
    df = pd.read_csv(tsv_file, sep='\t')
    original_row_count = len(df)

    print("Assigning orientation to hits...")
    df = add_orientation_column(df)

    print("Finding overlapping locus hits & keeping longest subject sequence per group...")
    df = by_subject_seq(df)
    final_count = len(df)

    print("Writing filtered hits to .tsv file...")
    df.to_csv(tsv_file.replace('_hits.tsv', '_distinct.tsv'), sep='\t', index=False)

    log_messages = [
        f"Original BLAST mappings: {original_row_count}",
        f"Count after overlap grouping & keeping largest sequence: {final_count}"
    ]

    create_log_file(tsv_file.replace('_hits.tsv', '_distinct.log'), log_messages)

    print("Processing completed. Check log file for details.")

if __name__ == "__main__":
    main()
