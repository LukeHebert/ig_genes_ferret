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
-query heavy/2_blast/J-3pRACE_HeavyChain_filtered.fasta \
-task megablast \
-db heavy/2_blast/IMGT_heavy_exons_db/IMGT_heavy_exons_fixed.fasta \
-out heavy/2_blast/J-3pRACE_HeavyChain_blastCs.tsv \
-outfmt "6 qseqid qlen sseqid slen qstart qend sstart send qseq bitscore evalue length" \
-num_threads 40
