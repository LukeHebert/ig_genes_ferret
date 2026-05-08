#!/usr/bin/env python3
"""
This script takes a ".assembled.fastq" file from merged BCR mRNA Illumina
sequencing as input and calls IgBLAST to map the V, D, and J genes within
sequences to their germline sequences.

Example use:
python identify_genes.py in.assembled.fastq ferret --igblast_path /path/to/igblastn --igblast_data_dir /path/to/igblast_data
"""

__author__ = "Luke S Hebert"
__version__ = "0.1.2"
__license__ = "MIT"

import argparse
import datetime
import os
import shlex
import subprocess
import sys
from Bio import SeqIO


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


def validate_directory(path, description):
    """Raise a clear error if a required directory does not exist."""
    if not os.path.isdir(path):
        raise NotADirectoryError(f"{description} was not found: {path}")
    return os.path.abspath(path)


def build_igblast_paths(igblast_data_dir, species):
    """Build the species-specific IgBLAST data paths."""
    data_dir = validate_directory(igblast_data_dir, "IgBLAST data directory")
    internal_dir = validate_directory(
        os.path.join(data_dir, "internal_data"),
        "IgBLAST internal_data directory",
    )
    optional_dir = validate_directory(
        os.path.join(data_dir, "optional_file"),
        "IgBLAST optional_file directory",
    )
    species_dir = validate_directory(
        os.path.join(internal_dir, species),
        f"IgBLAST species directory for '{species}'",
    )
    return {
        "v": validate_file(os.path.join(species_dir, f"{species}_V"), f"IgBLAST V database for '{species}'"),
        "j": validate_file(os.path.join(species_dir, f"{species}_J"), f"IgBLAST J database for '{species}'"),
        "d": validate_file(os.path.join(species_dir, f"{species}_D"), f"IgBLAST D database for '{species}'"),
        "c": validate_file(os.path.join(species_dir, f"{species}_C"), f"IgBLAST C database for '{species}'"),
        "aux": validate_file(os.path.join(optional_dir, f"{species}_gl.aux"), f"IgBLAST auxiliary file for '{species}'"),
    }


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
    """ Main entry point of the app """
    reads_path = validate_file(args.reads, "Input FASTQ file")
    igblast_path = validate_executable(args.igblast_path, "IgBLAST executable")
    igblast_data_paths = build_igblast_paths(args.igblast_data_dir, args.sp)

    now = datetime.datetime.now()
    directory = os.path.dirname(reads_path)
    log_name = f'log_mapping_{now.strftime("%Y-%m-%d_%H-%M-%S")}.txt'
    log_path = os.path.join(directory, log_name)

    args_str = call_igblast(reads_path, args.sp, log_path, igblast_path, igblast_data_paths)

    elapsed_str = time_passed(now)
    with open(log_path, 'a') as log:
        log.write(f'\n\nCOMMAND(s):\n{args_str}\n\n'
            f'MAPPING-ANNOTATING TIME\n(HR:MIN:SEC):\n{elapsed_str}')
    

def call_igblast(infile, species, log_file_name, igblast_path, igblast_data_paths):
    """Calls IgBLAST for mapping and annotating gene elements from BCR mRNA sequences."""
    print('\tConverting .fastq to .fasta...')
    infasta = infile.replace('.fastq','.fasta')
    sequences = SeqIO.parse(infile, "fastq")
    SeqIO.write(sequences, infasta, "fasta")

    outnametsv = infasta.replace('assembled.fasta','IgBLAST.tsv')
    command = [
        igblast_path,
        "-germline_db_V",
        igblast_data_paths["v"],
        "-germline_db_J",
        igblast_data_paths["j"],
        "-germline_db_D",
        igblast_data_paths["d"],
        "-c_region_db",
        igblast_data_paths["c"],
        "-auxiliary_data",
        igblast_data_paths["aux"],
        "-organism",
        species,
        "-query",
        infasta,
        "-outfmt",
        "19",
        "-out",
        outnametsv,
    ]
    print('\tCalling IgBLAST...')
    return run_and_log(command, log_file_name)


def time_passed(start_time):
    """ Makes a human-readable string of ellapsed time from start_time to
    when this function is called"""
    elapsed_time = datetime.datetime.now() - start_time
    elapsed_hr = int(elapsed_time.total_seconds() // 3600)
    elapsed_min = int((elapsed_time.total_seconds() % 3600) // 60)
    elapsed_sec = int(elapsed_time.total_seconds() % 60)
    return f"{elapsed_hr:02}:{elapsed_min:02}:{elapsed_sec:02}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=('Calls IgBLAST with preset parameters to annotate each read\'s genetic elements.'))
    parser.add_argument('reads', metavar='input_fastq', type=str,
        help='Input (assembled aka merged) .fastq file name.')
    parser.add_argument('sp', metavar='species', type=str,
        help='Species of sample donor e.g. "human", "mouse", or "ferret".')
    parser.add_argument(
        '--igblast_path',
        required=True,
        help='Path to the IgBLAST executable (igblastn).',
    )
    parser.add_argument(
        '--igblast_data_dir',
        required=True,
        help='Base directory containing IgBLAST internal_data and optional_file directories.',
    )

    args = parser.parse_args()
    
    try:
        main(args)
    except (FileNotFoundError, PermissionError, NotADirectoryError) as error:
        print(error, file=sys.stderr)
        sys.exit(1)
