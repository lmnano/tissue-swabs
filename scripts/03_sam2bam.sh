#!/bin/bash

samtools view -bS /path/to/bowtie/sam/file.sam > /path/to/bam/output.bam
samtools sort /path/to/bam/output.bam > /path/to/sorted/bam/output_sorted.bam
samtools index /path/to/sorted/bam/output_sorted.bam /path/to/bai/index/output.bai
