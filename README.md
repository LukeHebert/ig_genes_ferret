# ig_genes_ferret
A detailed workflow of how I (and colleagues) predicted _Mustela mutorius furo_ (the domestic ferret) immunoglobulin (Ig) genes using the ferret genome assembly ASM1176430v1.1 and then mapped ferret B cell receptor sequencing data to those predictions.

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
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/) (version 4.6)
- [PEAR](https://github.com/tseemann/PEAR) (version 0.9.11)
- [SignalP](https://services.healthtech.dtu.dk/services/SignalP-5.0/) (version 5.0b)
