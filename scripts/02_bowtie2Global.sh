#!/bin/bash

bowtie2 --very-sensitive \
        -x /path/to/bt2/index \
        --threads 8 \
        --met-file bt2mterics.txt \
        --al-gz /path/to/aligned/output/reads/file.fq.gz \
        --un-gz /path/to/unaligned/reads/file.fq.gz \
        -U /path/to/catalog/reads/fasta/file.fa.gz \
        -f \
        -t \
        -S /path/to/output/sam/file.sam;
