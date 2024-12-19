# ig_genes_ferret
A detailed workflow of how I (and colleagues) predicted _Mustela mutorius furo_ (the domestic ferret) immunoglobulin (Ig) genes using the ferret genome assembly ASM1176430v1.1 and then mapped ferret B cell receptor sequencing data to those predictions.

# Table of contents
- `variable_gene_prediction/`: A directory of all scripts used in the prediction of ferret variable region Ig genes (i.e. IGHV, IGHD, IGHJ, IGKV, IGKJ, IGLV, IGLJ gene segments) and a description of how they were used. These scripts are generalizable to the discovery of Ig genes in other species, though some customization is likely necessary.
- 
- `variable_gene_transcript_analysis/` A directory of all scripts used to analyze B cell receptor variable domain-encoding mRNA (cDNA) transcription data from ferrets, and a description of how they were used. This analysis uses the gene sequences obtained with scripts from the variable_gene_prediction directory.
- 
- `constant_gene_transcript_analysis/` A directory of all scripts used to analyze B cell receptor constant domain-encoding mRNA (cDNA) transcription data from ferrets, and a description of how they were used.
