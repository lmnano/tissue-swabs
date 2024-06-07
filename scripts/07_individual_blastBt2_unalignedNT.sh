#!/bin/bash

FILESDIR=/path/to/files/dir
OUTDIR=/path/to/out/dir

#if running all the files from the directory
FILES=$(ls $FILESDIR | cat)

cd $FILESDIR

#if running only a subset of files
#FILES2=$(ls PC639_s.fasta PC639_t.fasta PC641_s.fasta PC641_t.fasta PC680_s.fasta | cat)


for f in $FILES2;
	do

		OUTF=$(echo "$f" | cut -d '.' -f 1)

		blastn -db /path/to/nt/database/nt \
			-query "$f" \
			-out $OUTDIR/$OUTF.limit.asn1 \
			-num_threads 10 \
			-evalue 0.001 \
			-perc_identity 50 \
			-outfmt 11;

		blast_formatter -archive $OUTDIR/$OUTF.limit.asn1 -outfmt 0 -out $OUTDIR/$OUTF.limit.outfmt0.txt
	done;

