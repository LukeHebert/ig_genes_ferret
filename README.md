# ig_genes_ferret
A detailed workflow of how I (and colleagues) predicted _Mustela mutorius furo_ (the domestic ferret) immunoglobulin (Ig) genes using the ferret genome assembly ASM1176430v1.1 and then mapped ferret B cell receptor sequencing data to those predictions.

## Important note
This repository is best suited for reproducing the results associated with the author's ferret immunoglobulin gene publication:
[Hebert et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12908631/#_ad93_).

If you are trying to annotate the immunoglobulin genes of a new genome assembly, especially a mammalian assembly, I recommend first trying one or both of these dedicated tools and their associated publications:
- [IgDetective](https://github.com/Immunotools/IgDetective)
- [Digger](https://williamdlees.github.io/digger/_build/html/index.html)

If you are only looking for the author's ferret germline reference sequences and their associated metadata, those are best obtained from `DataS4_ferret_ig_gene_annotations` in the supplementary materials of the same publication:
[Hebert et al. supplementary materials](https://pmc.ncbi.nlm.nih.gov/articles/PMC12908631/#_ad93_:~:text=in%20this%20article.-,Supplementary%20Materials,-DataS3_ferret_J_alleles).

# Table of contents
- `variable_gene_prediction/`: A directory containing all scripts used in the prediction of ferret variable region Ig genes (i.e. IGHV, IGHD, IGHJ, IGKV, IGKJ, IGLV, IGLJ gene segments) and a description of how they were used. These scripts are generalizable to the discovery of Ig genes in other species, though some customization is likely necessary.
  
- `variable_gene_transcript_analysis/` A directory containing all scripts used to analyze B cell receptor variable domain-encoding mRNA (cDNA) transcription data from ferrets, and a description of how they were used. This analysis uses the gene sequences obtained with scripts from the "variable_gene_prediction/" directory.
  
- `constant_gene_transcript_analysis/` A directory containing all scripts used to analyze B cell receptor constant domain-encoding mRNA (cDNA) transcription data from ferrets, and a description of how they were used.

# Dependencies

### Python packages
All python dependencies are listed in the `environment.yml` file.

### Other tools
Non-python package bioinformatics tools called by the scripts in this repository include:
- [BLAST+](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html#blast-executables) suite (version 2.15.0)
- [Digger](https://williamdlees.github.io/digger/_build/html/index.html) (version 0.7.3)
- [IgBLAST](https://ncbi.github.io/igblast/) (version 1.21.0)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (version 4.6)
- [PEAR](https://github.com/tseemann/PEAR) (version 0.9.11)
- [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/) (version 5.0b)

The scripts in this repository now expect you to provide the paths to these external tools on the command line, rather than storing machine-specific tool paths in the code.
