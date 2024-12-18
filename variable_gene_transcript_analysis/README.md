# Most important scripts used in the analysis of variable domain-encoding gene expression in the ferret

## A. Filtering and mapping
1. *Quality control filtering and merging:* The script `trim_merge.py` was used on all B cell receptor (variable gene-containing) Illumina paired-end read (FASTQ.GZ) files. This script calls a local instantiation of [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [PEAR](https://github.com/tseemann/PEAR).
2. *Mapping reads to reference Ig genes:* The script `identify_genes.py` was used to map & align cDNA reads to a reference set of ferret genes by calling a local instance of [IgBLAST](https://ncbi.github.io/igblast/cook/How-to-set-up.html). The reference gene dataset was formatted for IgBLAST with the help of [IMGT's V-QUEST](https://www.imgt.org/IMGT_vquest/input) (especially for the creation of `.ndm.imgt` and `_gl.aux` files).

## B. Assessing gene usage
1. The outputs from step A2 were passed to the scripts `heatmap_tsv.py` and then `heatmap_plot.py` to generate gene usage statistics and plots, respectively, from all variable gene amplicons.

## C. Analyzing CDR3 sequences
1. The three scripts `cdr3_aa_distro.py`, `cdr3_hydro.py`, and `cdr3_lengths.py` were used in no particular order to generate descriptive (*in-silico* translated) CDR3 amino acid statistics and plots. As the names imply, `cdr3_aa_distro.py` analyzes amino acid usage, `cdr3_hydro.py` analyzes hydrophobicity of CDR3s, and `cdr3_lengths.py` analyzes CDR3 lengths.
