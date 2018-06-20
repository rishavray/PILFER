#!/usr/bin/python

import argparse
import os
import gzip
import sys

#Defining the argparse object and its attributes
parser = argparse.ArgumentParser(description="Collapses the reads by merging the same reads with counts as the fasta header")
requiredArgument = parser.add_argument_group('Required arguments')
requiredArgument.add_argument("-i", metavar="<input filename>", help="Input fastq file. - if input is stdin",required=True)
requiredArgument.add_argument("-o", metavar="<output filename>", help="Output fastq file. - if output is stdout",required=True)
requiredArgument.add_argument("--format", metavar="<format>", choices=['fastq', 'fasta'], help="Input format. Must be 'fastq' or 'fasta'",required=True)
args = parser.parse_args()

#Setting attributes and Sanity check

infile = args.i
outfile = args.o
in_format = args.format

if infile == "-":
	infile_pt = sys.stdin
	sys.stderr.write("Reading from stdin\n")
else:
	if not os.path.isfile(infile):
		sys.stderr.write("Not a valid input file\n")
		exit()
	else:
		sys.stderr.write("Reading from file " + infile + "\n")
		if infile.endswith(".gz"):
			infile_pt = gzip.open(infile,"rb")
		else:
			infile_pt = open(infile,"rb")

if outfile == "-":
	outfile_pt = sys.stdout
else:
	if not outfile.endswith(".gz"):
		outfile += ".gz"
	outfile_pt = gzip.open(outfile,"wb")

#Collapsing reads

seq_dict = {}
sys.stderr.write("Processing reads\n")
original_read_count = 0
while True:
	line = infile_pt.readline().strip()
	if not line:
		break
	seq = infile_pt.readline().strip()

	
	if in_format == "fastq":
		
		if seq in seq_dict:
			seq_dict[seq] += 1
		else:
			seq_dict[seq] = 1
		infile_pt.readline().strip()
		infile_pt.readline().strip()
	else:

		if seq in seq_dict:
			seq_dict[seq] += 1
		else:
			seq_dict[seq] = 1

	original_read_count += 1

#printing collapsed reads
i = 0
for key in seq_dict:
	seq_id = 'seq' + str(i).zfill(6)
	outfile_pt.write(">" + str(seq_dict[key]) + ":" + seq_id + "\n" + key + "\n")
	i += 1

outfile_pt.close()
infile_pt.close()
sys.stderr.write("Original number of reads: " + str(original_read_count) + "\n")
sys.stderr.write("collapsed number of reads: " + str(len(seq_dict)) + "\n")