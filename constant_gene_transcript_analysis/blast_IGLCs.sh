#!/bin/bash

BLASTN="$1"

if [ -z "$BLASTN" ]; then
    echo "Usage: $0 /path/to/blastn" >&2
    exit 1
fi

if [ ! -x "$BLASTN" ]; then
    echo "BLASTN executable is missing or not executable: $BLASTN" >&2
    exit 1
fi

"$BLASTN" \
-query lambda/2_blast/F2_IGLC_merged.fasta \
-task megablast \
-db lambda/2_blast/IMGT_lambda_Cexons_db/IMGT_lambda_Cexons_fixed.fasta \
-out lambda/2_blast/F2_IGLC_merged.tsv \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq bitscore evalue length" \
-num_threads 40
