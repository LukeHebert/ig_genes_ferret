/stor/work/Georgiou/Sharing_Folder/BLASTplus_2.15.0/ncbi-blast-2.15.0+/bin/blastn \
-query kappa/2_blast/F2_IGKC_merged.fasta \
-task megablast \
-db kappa/2_blast/IMGT_kappa_Cexons_db/IMGT_kappa_Cexons_fixed.fasta \
-out kappa/2_blast/F2_IGKC_merged.tsv \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq bitscore evalue length" \
-num_threads 40