#!/bin/bash

makeblastdb -in /path/to/fasta/file \
	-dbtype nucl \
	-out /path/to/out/folder \
	-parse_seqids
