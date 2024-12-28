# Analyzing transcript data of B cell receptor constant genes

## A. Quality control filtering
1. <ins>**Heavy chain data:**</ins> After amplifying B cell receptor transcripts with a J gene forward primer pool and 3'RACE, the amplicon was sequenced using Oxford Nanopore (Sequencing Kit V14 and the NEB Next Companion Module + Dorado SUP basecaller). Both 5' and 3' sequencing ends were trimmed with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and any sequences with an average Phred score below 20 were removed. This was all done by calling the script on the merged ONT data: `qc_nanopore.py`.
   
2. <ins>**Light chain data:**</ins> Light chain C gene transcript data was sequenced with Illumina MiSeq 2 x 300 (paired-end), which was filtered for quality control and merged with the same scripts used for variable gene data processing, `trim_merge.py`. This script calls local instances of [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) and [PEAR](https://github.com/tseemann/PEAR).

## B. Assessing gene usage
1. <ins>**BLAST querying C gene exons:**</ins> First, a custom, local [BLAST database was created](https://www.ncbi.nlm.nih.gov/books/NBK569841/) using previously-reported [ferret C gene data from IMGT](https://www.imgt.org/genedb/). Exon combinations were used such that each C gene exon could be mapped to individually or combinations could be mapped to, depending on the length and identity of a given read. This subverts the primary challenge of using BLAST for transcript analysis (namely, the inability to account for intronic sequences in a reference database). The BLAST commands used for this analysis are recorded in `blast_IGHCs.sh`, `blast_IGKCs`, and `blast_IGLCs`.
   
2. <ins>**Filtering & plotting gene usage:**</ins> Next, `plot_lengths_seqids.py` was called; in this script, BLAST hits were filtered so that only the top hit per (transcript) query was kept, based on bitscore. Last, the remaining hits were plotted for their mappings (both generally as C gene usage and specifically as exon combination identity) and alignment lengths.
