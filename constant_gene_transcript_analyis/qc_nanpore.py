#!/usr/bin/env python3

import argparse
import os
import subprocess
import datetime
from Bio import SeqIO

def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter FASTQ reads by trimming low-quality bases at both ends using cutadapt, "
                    "then remove reads whose average quality falls below a specified threshold."
    )
    parser.add_argument("fastq", help="Path to the input FASTQ file.")
    parser.add_argument(
        "--trim_quality",
        type=int,
        required=True,
        help="Quality cutoff used by cutadapt for trimming at both ends."
    )
    parser.add_argument(
        "--avg_quality",
        type=int,
        required=True,
        help="Minimum average quality threshold to retain reads after trimming."
    )
    return parser.parse_args()


def run_cutadapt(input_fastq, trim_quality, output_trimmed_fastq):
    """Run cutadapt to trim reads at both ends based on the specified quality cutoff."""
    cutadapt_path = "/usr/local/bin/cutadapt"
    cmd = [
        cutadapt_path,
        "-q", str(trim_quality),    # One value applies trimming to both 5' and 3' ends
        "-j", "40",
        "-o", output_trimmed_fastq,
        input_fastq
    ]
    subprocess.run(cmd, check=True)
    return cmd


def filter_by_average_quality(input_fastq, output_fastq, avg_quality_cutoff):
    retained_count = 0
    with open(output_fastq, "w") as out_handle:
        for record in SeqIO.parse(input_fastq, "fastq"):
            # Skip empty reads
            if len(record) == 0:
                continue

            avg_qual = sum(record.letter_annotations["phred_quality"]) / len(record)
            if avg_qual >= avg_quality_cutoff:
                SeqIO.write(record, out_handle, "fastq")
                retained_count += 1
    return retained_count



def count_reads_in_fastq(fastq_path):
    """Count the number of reads in a FASTQ file."""
    return sum(1 for _ in SeqIO.parse(fastq_path, "fastq"))


def write_log(log_file, cmd, input_fastq, trimmed_fastq, filtered_fastq, trim_quality, avg_quality, original_count, trimmed_count, filtered_count):
    """Write a detailed log file."""
    with open(log_file, "w") as lf:
        lf.write("Command run:\n")
        lf.write(" ".join(cmd) + "\n\n")
        lf.write(f"Date/Time: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        lf.write(f"Input FASTQ: {input_fastq}\n")
        lf.write(f"Trimmed FASTQ: {trimmed_fastq}\n")
        lf.write(f"Final Filtered FASTQ: {filtered_fastq}\n")
        lf.write(f"Trim quality cutoff: {trim_quality}\n")
        lf.write(f"Average quality cutoff: {avg_quality}\n")
        lf.write(f"Number of reads (original): {original_count}\n")
        lf.write(f"Number of reads (after trimming): {trimmed_count}\n")
        lf.write(f"Number of reads (after filtering): {filtered_count}\n")
        lf.write("Cutadapt arguments used:\n")
        lf.write(" ".join(cmd) + "\n")


def main():
    args = parse_args()

    input_fastq = args.fastq
    trim_quality = args.trim_quality
    avg_quality = args.avg_quality

    # Validate input
    if not os.path.isfile(input_fastq):
        raise FileNotFoundError(f"Input FASTQ file '{input_fastq}' does not exist.")

    # Paths and naming
    output_dir = os.path.dirname(os.path.abspath(input_fastq))
    base_name = os.path.basename(input_fastq)
    root_name, ext = os.path.splitext(base_name)

    # Intermediate and final output files
    trimmed_fastq = os.path.join(output_dir, f"{root_name}_trimmed{ext}")
    filtered_fastq = os.path.join(output_dir, f"{root_name}_filtered{ext}")
    now = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = os.path.join(output_dir, f"{root_name}_filter_log_{now}.txt")

    # Count original reads
    original_count = count_reads_in_fastq(input_fastq)

    # Run cutadapt
    cmd = run_cutadapt(input_fastq, trim_quality, trimmed_fastq)

    # Count reads after trimming
    trimmed_count = count_reads_in_fastq(trimmed_fastq)

    # Filter by average quality
    filtered_count = filter_by_average_quality(trimmed_fastq, filtered_fastq, avg_quality)

    # Write log file
    write_log(log_file, cmd, input_fastq, trimmed_fastq, filtered_fastq, trim_quality, avg_quality,
              original_count, trimmed_count, filtered_count)


if __name__ == "__main__":
    main()
