# Combine distinct BLAST hits + annotated adjacent sequences, while only keeping one row of column headers
cat genes_v/IGxV_seqs.tsv > all_distinct_seqs.tsv
tail -n +2 genes_d/IGHD_seqs.tsv >> all_distinct_seqs.tsv
tail -n +2 genes_j/IGxJ_seqs.tsv >> all_distinct_seqs.tsv