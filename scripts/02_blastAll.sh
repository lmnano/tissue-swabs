#!/bin/bash

blastn -query /path/to/query/fasta \
	-db /path/to/blast/db \
	-outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sscinames staxids stitle sseq" \
	-num_threads 4 \
	-max_target_seqs 5

