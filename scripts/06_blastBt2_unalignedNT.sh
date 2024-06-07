#!/bin/bash

blastn -db  /path/to/nt/database \
	-query /path/to/unaligned/bowtie2/results \
	-out /path/to/output/file.asn1 \
	-num_threads 10 \
	-evalue 0.001 \
	-perc_identity 50 \
	-outfmt 11;

blast_formatter -archive /path/to/output/file.asn1 -outfmt 0 -out /path/to/outfmt0/output/file.outfmt0.txt

