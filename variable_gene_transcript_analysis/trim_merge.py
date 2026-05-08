#!/usr/bin/env python3
"""
Takes one forward and one reverse read file from Illumina platform paired-end
sequencing of B cell receptors and:
1) generates fastQC .html files to manually evaluate quality metrics of sequencing data
2) trims library prep adapter sequences & trims poor quality 3' ends of reads
3) merges the paired forward & reverse reads

A log file will be created in your current working directory.

Example use:
python trim_merge.py r1.fastq.gz r2.fastq.gz --fastqc_path /path/to/fastqc --cutadapt_path /path/to/cutadapt --pear_path /path/to/pear
"""

__author__ = "Luke S Hebert"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import datetime
import os
import shlex
import subprocess
import sys


def validate_file(path, description):
    """Raise a clear error if a required file does not exist."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"{description} was not found: {path}")
    return os.path.abspath(path)


def validate_executable(path, description):
    """Raise a clear error if a required executable path is invalid."""
    executable_path = validate_file(path, description)
    if not os.access(executable_path, os.X_OK):
        raise PermissionError(f"{description} is not executable: {executable_path}")
    return executable_path


def run_and_log(command, log_pathway):
    """Run a command, append stdout and stderr to the log, and return the command string."""
    command_str = shlex.join(command)
    with open(log_pathway, "a") as log_handle:
        log_handle.write(f"\n$ {command_str}\n")
        subprocess.run(
            command,
            check=True,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
        )
    return command_str


def main(args):
    """ Main entry point of the module """
    r1_path = validate_file(args.R1, "Forward reads file")
    r2_path = validate_file(args.R2, "Reverse reads file")
    fastqc_path = validate_executable(args.fastqc_path, "FastQC executable")
    cutadapt_path = validate_executable(args.cutadapt_path, "Cutadapt executable")
    pear_path = validate_executable(args.pear_path, "PEAR executable")

    # Create log file with current date and time in its name
    now = datetime.datetime.now()
    directory = os.path.dirname(r1_path)
    log_name = f'log_QC-trimming-merging_{now.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
    log_path = os.path.join(directory, log_name)

    # Run FastQC on initial input files
    qc_args = run_fastqc([r1_path, r2_path], directory, log_path, fastqc_path)
    
    if not args.notrim:
        # Trim adapters & poor quality read ends
        trimmed_r1, trimmed_r2, trim_args = call_cutadapt(
            r1_path,
            r2_path,
            directory,
            log_path,
            cutadapt_path,
        )
    else:
        # Or use original files if user specifies to skip trimming
        trimmed_r1, trimmed_r2 = r1_path, r2_path
        trim_args = "Trimming skipped due to --notrim flag"
    
    # Call "PEAR" Paired End reAd mergeR to consolidate forward R1 & reverse R2 reads
    merge_args = call_pear(trimmed_r1, trimmed_r2, log_path, pear_path)

    # Run FastQC on the output assembled file
    out_base = trimmed_r1.replace('_trimmed.fastq.gz', '').replace('_R1','')
    assembled_file = f"{out_base}.assembled.fastq"
    run_fastqc([assembled_file], directory, log_path, fastqc_path)
    
    # Append the tool-calling command and the time it took to log file
    elapsed_str = time_passed(now)
    with open(log_path, 'a') as log:
        log.write(f'\n\nQUALITY CHECK COMMAND:\n{qc_args}'
                  f'\n\nADAPTER & PHRED TRIMMING COMMAND:\n{trim_args}'
                  f'\n\nMERGE COMMAND:\n{merge_args}'
                  f'\n\nTOTAL SCRIPT RUN TIME\n(HR:MIN:SEC):\n{elapsed_str}')


def run_fastqc(files, output_dir, log_pathway, fastqc_path):
    """ Runs FastQC on a list of files and outputs the reports to the specified directory """
    command = [fastqc_path, "-o", output_dir] + files
    return run_and_log(command, log_pathway)


def call_cutadapt(r1, r2, output_dir, log_pathway, cutadapt_path):
    """ Trims adapters and low-quality ends from reads using cutadapt """
    #These are the TruSeq standard adapter sequences
    adapter1 = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCA'
    adapter2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT'
    quality_cutoff = '20'  # Quality score cutoff for 3' ends
    trimmed_r1 = os.path.join(output_dir, os.path.basename(r1).replace('.fastq.gz', '_trimmed.fastq.gz'))
    trimmed_r2 = os.path.join(output_dir, os.path.basename(r2).replace('.fastq.gz', '_trimmed.fastq.gz'))
    # Note: -m 50 (keep minimum length cutoff) is needed to prevent downstream PEAR merging errors with sequence records of length 0
    # Note: the --pair-filter=both is needed to prevent uneven output file lengths which also causes a downstream PEAR error
    cmd = [
        cutadapt_path,
        "-q",
        quality_cutoff,
        "-a",
        adapter1,
        "-A",
        adapter2,
        "-o",
        trimmed_r1,
        "-p",
        trimmed_r2,
        r1,
        r2,
        "-m",
        "50",
        "-j",
        "0",
    ]
    command_str = run_and_log(cmd, log_pathway)
    run_and_log(cmd, log_pathway)
    return trimmed_r1, trimmed_r2, command_str


def call_pear(r1, r2, log_pathway, pear_path):
    """ Calls PEAR read merger to find forward-reverse read couples based on
    overlapping sequence ends """
    
    out = r1.replace('_trimmed.fastq.gz', '').replace('_R1','')
    command = [
        pear_path,
        "-f",
        r1,
        "-r",
        r2,
        "-o",
        out,
        "-v",
        "10",
        "-m",
        "700",
        "-n",
        "50",
        "-u",
        "1",
        "-j",
        "20",
    ]
    return run_and_log(command, log_pathway)


def time_passed(start_time):
    """ Makes a human-readable string of elapsed time from start_time to 
    when this function is called"""
    elapsed_time = datetime.datetime.now() - start_time
    elapsed_hr = int(elapsed_time.total_seconds() // 3600)
    elapsed_min = int((elapsed_time.total_seconds() % 3600) // 60)
    elapsed_sec = int(elapsed_time.total_seconds() % 60)
    return f"{elapsed_hr:02}:{elapsed_min:02}:{elapsed_sec:02}"


if __name__ == "__main__":
    """ This is executed when run from the command line """
    
    # Get the arguments from terminal
    parser = argparse.ArgumentParser(
        description='Generates .html files for quality assessment, trims adapters and low quality 3 prime read ends, and merges paired forward & reverse reads.')
    parser.add_argument('R1', metavar='forward_reads', type=str, 
        help='R1 fastq.gz file')
    parser.add_argument('R2', metavar='reverse_reads', type=str, 
        help='R2 fastq.gz file')
    parser.add_argument('--notrim', action='store_true', 
                        help='If set, skips the adapter & quality score trimming step before merging paired reads.')
    parser.add_argument('--fastqc_path', required=True, help='Path to the FastQC executable.')
    parser.add_argument('--cutadapt_path', required=True, help='Path to the Cutadapt executable.')
    parser.add_argument('--pear_path', required=True, help='Path to the PEAR executable.')
    args = parser.parse_args()
    
    # Pass the arguments to the main pairing function
    try:
        main(args)
    except (FileNotFoundError, PermissionError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
