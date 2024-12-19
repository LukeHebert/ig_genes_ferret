/stor/work/Georgiou/Sharing_Folder/BLASTplus_2.15.0/ncbi-blast-2.15.0+/bin/blastn \
-query lambda/2_blast/F2_IGLC_merged.fasta \
-task megablast \
-db lambda/2_blast/IMGT_lambda_Cexons_db/IMGT_lambda_Cexons_fixed.fasta \
-out lambda/2_blast/F2_IGLC_merged.tsv \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq bitscore evalue length" \
-num_threads 40