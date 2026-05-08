# Analysis of variable domain-encoding Ig gene expression in the ferret

## A. Filtering and mapping
1. <ins>**Quality control filtering and merging:**</ins> The script `trim_merge.py` was used on all B cell receptor (variable gene-containing) Illumina paired-end read (FASTQ.GZ) files. You now provide the paths to [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), and [PEAR](https://github.com/tseemann/PEAR) on the command line, for example:

`python trim_merge.py reads_R1.fastq.gz reads_R2.fastq.gz --fastqc_path /path/to/fastqc --cutadapt_path /path/to/cutadapt --pear_path /path/to/pear`
  
2. <ins>**Mapping reads to reference Ig genes:**</ins> The script `identify_genes.py` was used to map & align cDNA reads to a reference set of ferret genes by calling [IgBLAST](https://ncbi.github.io/igblast/cook/How-to-set-up.html). You now provide both the `igblastn` executable path and the base IgBLAST data directory on the command line, for example:

`python identify_genes.py sample.assembled.fastq ferret --igblast_path /path/to/igblastn --igblast_data_dir /path/to/igblast_data`

The reference gene dataset was formatted for IgBLAST with the help of [IMGT's V-QUEST](https://www.imgt.org/IMGT_vquest/input) (especially for the creation of `.ndm.imgt` and `_gl.aux` files).

## B. Assessing gene usage
1. The outputs from step A2 were passed to the scripts `heatmap_tsv.py` and then `heatmap_plot.py` to generate gene usage statistics and plots, respectively, from all variable gene amplicons.

## C. Analyzing CDR3 sequences
1. The three scripts `cdr3_aa_distro.py`, `cdr3_hydro.py`, and `cdr3_lengths.py` were used in no particular order to generate descriptive (*in-silico* translated) CDR3 amino acid statistics and plots. As the names imply, `cdr3_aa_distro.py` analyzes amino acid usage, `cdr3_hydro.py` analyzes hydrophobicity of CDR3s, and `cdr3_lengths.py` analyzes CDR3 lengths.
