#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import datetime
import statistics

def plot_length_hist(df, base_filename, plot_color):
    """
    Plot the distribution of the 'length' column values as a histogram.
    """
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df, x="length", kde=False, color=plot_color)
    plt.title("Distribution of Length")
    plt.xlabel("Length")
    plt.ylabel("Count")
    output_filename = f"{base_filename}_length_hist.png"
    plt.savefig(output_filename, dpi=800)
    plt.close()

def plot_sseqid_bar(df, base_filename, plot_color):
    """
    Plot the 'sseqid' column values as a bar plot (countplot).
    Each unique sseqid is one bar.
    """
    plt.figure(figsize=(10, 6))
    order = df["sseqid"].value_counts().index
    sns.countplot(y="sseqid", data=df, order=order, color=plot_color)
    plt.title("Count of Each sseqid")
    plt.xlabel("Count")
    plt.ylabel("sseqid")
    output_filename = f"{base_filename}_sseqid_bar.png"
    plt.savefig(output_filename, dpi=800, bbox_inches='tight')
    plt.close()
    return order, df["sseqid"].value_counts()

def plot_sseqid_bar_parsed(df, base_filename, plot_color):
    """
    Parse the 'sseqid' column by splitting on '_' and taking the first element,
    then plot as a bar plot (countplot).
    """
    df["sseqid_parsed"] = df["sseqid"].str.split("_").str[0]
    plt.figure(figsize=(10, 6))
    order = df["sseqid_parsed"].value_counts().index
    sns.countplot(y="sseqid_parsed", data=df, order=order, color=plot_color)
    plt.title("Count of Parsed sseqid (First Segment Before '_')")
    plt.xlabel("Count")
    plt.ylabel("Parsed sseqid")
    output_filename = f"{base_filename}_sseqid_parsed_bar.png"
    plt.savefig(output_filename, dpi=800, bbox_inches='tight')
    plt.close()
    return order, df["sseqid_parsed"].value_counts()

def main():
    parser = argparse.ArgumentParser(description="Plot histograms and barplots from a TSV file with BLASTn filtered results.")
    parser.add_argument("tsv_file", help="The input TSV file (with headers including 'sseqid', 'qseqid', 'bitscore', and 'length').")
    parser.add_argument("--color", help="Optional color/hexcode for plots.", default="#4c72b0")
    args = parser.parse_args()

    # Load data
    df = pd.read_csv(args.tsv_file, sep="\t")
    
    # For any group of rows that share a qseqid, keep only the row with the top bitscore
    df = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]

    # Derive a base filename for output
    base_filename = os.path.splitext(args.tsv_file)[0]

    # Save filtered dataframe
    df.to_csv(f"{base_filename}_topBitPerQuery.tsv", sep='\t')

    # Create plots
    plot_length_hist(df, base_filename, args.color)
    sseqid_order, sseqid_counts = plot_sseqid_bar(df, base_filename, args.color)
    sseqid_parsed_order, sseqid_parsed_counts = plot_sseqid_bar_parsed(df, base_filename, args.color)

    # Logging
    # Create a datetime-stamped log file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{base_filename}_{timestamp}.log"

    # Collect stats
    length_mean = df["length"].mean()
    length_std = df["length"].std()
    length_min = df["length"].min()
    length_max = df["length"].max()
    length_count = df["length"].count()

    # Write log
    with open(log_filename, "w") as log_f:
        # First line: user's terminal command
        log_f.write("Command used:\n")
        log_f.write(" ".join(sys.argv) + "\n\n")

        # Report how many rows contributed
        log_f.write(f"Number of rows contributing to the plots: {len(df)}\n\n")

        # A) Stats for the length plot
        log_f.write("Length column statistics:\n")
        log_f.write(f" Mean: {length_mean}\n")
        log_f.write(f" Std Dev: {length_std}\n")
        log_f.write(f" Min: {length_min}\n")
        log_f.write(f" Max: {length_max}\n")
        log_f.write(f" Count: {length_count}\n\n")

        # B) For each bar plot, report x tick labels and their corresponding heights
        # For the sseqid bar plot (y = sseqid, x = counts)
        log_f.write("sseqid Bar Plot Data (label: count):\n")
        for label in sseqid_order:
            log_f.write(f" {label}: {sseqid_counts[label]}\n")
        log_f.write("\n")

        # For the parsed sseqid bar plot
        log_f.write("sseqid_parsed Bar Plot Data (label: count):\n")
        for label in sseqid_parsed_order:
            log_f.write(f" {label}: {sseqid_parsed_counts[label]}\n")

if __name__ == "__main__":
    main()
