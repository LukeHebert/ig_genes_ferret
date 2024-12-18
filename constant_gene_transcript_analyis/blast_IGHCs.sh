/stor/work/Georgiou/Sharing_Folder/BLASTplus_2.15.0/ncbi-blast-2.15.0+/bin/blastn \
-query heavy/2_blast/J-3pRACE_HeavyChain_filtered.fasta \
-task megablast \
-db heavy/2_blast/IMGT_heavy_exons_db/IMGT_heavy_exons_fixed.fasta \
-out heavy/2_blast/J-3pRACE_HeavyChain_blastCs.tsv \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq bitscore evalue length" \
-num_threads 40