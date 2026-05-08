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
-query kappa/2_blast/F2_IGKC_merged.fasta \
-task megablast \
-db kappa/2_blast/IMGT_kappa_Cexons_db/IMGT_kappa_Cexons_fixed.fasta \
-out kappa/2_blast/F2_IGKC_merged.tsv \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq bitscore evalue length" \
-num_threads 40
