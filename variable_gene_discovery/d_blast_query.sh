#!/bin/bash

# Define variables
BLASTN="/stor/work/Georgiou/Sharing_Folder/BLASTplus_2.15.0/ncbi-blast-2.15.0+/bin/blastn"
QUERY_FILE="IGHD_all.fasta"
DB_NAME="../../ferret_blast_db/ferret_blastable"
OUTPUT_FILE="IGHD_hits.tsv"
EVALUE="1e2"
MAX_TARGET_SEQS="200"
DUST="no"
SOFT_MASKING="false"
PENALTY="-3"
REWARD="1"
GAP_OPEN="5"
GAP_EXTEND="2"
STRAND="both"
WORD_SIZE="7"
OUT_FORMAT="6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
HEADER="query acc.ver\tsubject acc.ver\t% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score"

# Run BLASTn
$BLASTN -query $QUERY_FILE -db $DB_NAME -outfmt "$OUT_FORMAT" -evalue $EVALUE -max_target_seqs $MAX_TARGET_SEQS -dust $DUST -soft_masking $SOFT_MASKING -penalty $PENALTY -reward $REWARD -gapopen $GAP_OPEN -gapextend $GAP_EXTEND -strand $STRAND -word_size $WORD_SIZE -out $OUTPUT_FILE

# Add headers to the output file
echo -e $HEADER | cat - $OUTPUT_FILE > temp.tsv && mv temp.tsv $OUTPUT_FILE
