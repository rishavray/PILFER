#!/usr/bin/python

import argparse
import os
import gzip
import sys
import uuid

#Function for writing reads
def writeRead(read,out_pointer,out_format):
	if out_format == "fastq":
		out_pointer.write(read[0] + "\n" + read[1] + "\n" + read[2] + "\n" + read[3] + "\n")
	else:
		out_pointer.write(read[0] + "\n" + read[1] + "\n")


#Defining the argparse object and its attributes
parser = argparse.ArgumentParser(description="Filters read with repeated motifs")
parser.add_argument("--min", default=3, metavar="<value>", type=int, help="Minimum size of the motif <value>")
parser.add_argument("--max", default=4, metavar="<value>", type=int, help="Maximum size of the motif <value>")
parser.add_argument("--threshold", default=0.75, metavar="<value>", type=float, help="Minimum size of the motif <value>")
requiredArgument = parser.add_argument_group('Required arguments')
requiredArgument.add_argument("-i", metavar="<input filename>", help="Input fastq file. - if input is stdin",required=True)
requiredArgument.add_argument("-o", metavar="<output filename>", help="Output fastq file. - if output is stdout",required=True)
requiredArgument.add_argument("--format", metavar="<format>", choices=['fastq', 'fasta'], help="Input format. Must be 'fastq' or 'fasta'",required=True)
args = parser.parse_args()

#setting attributes and sanity check
infile = args.i
outfile = args.o
in_format = args.format
threshold = args.threshold
min_motif = args.min
max_motif = args.max

if infile == "-":
	infile_pt = sys.stdin
	sys.stderr.write("Reading from stdin\n")
else:
	if not os.path.isfile(infile):
		sys.stderr.write("Not a valid input file\n")
		exit()
	else:
		sys.stderr.write("Reading from file" + infile + "\n")
		if infile.endswith(".gz"):
			infile_pt = gzip.open(infile,"rb")
		else:
			infile_pt = open(infile,"rb")

if outfile == "-":
	outfile_good_pt = sys.stdout
	if infile == "-":
		outfile_bad_pt = gzip.open(uuid.uuid4().get_hex()+"_bad.gz","wb")
	else:
		outfile_bad_pt = gzip.open(infile+"_bad.gz","wb")
else:
	if not outfile.endswith(".gz"):
		outfile_good_pt = gzip.open(outfile+"_good.gz","wb")
		outfile_bad_pt = gzip.open(outfile+"_bad.gz","wb")

dusty = 0
non_dusty = 0
sys.stderr.write("Processing reads\n")

while True:
	line = infile_pt.readline().strip()
	if not line:
		break
	seq = infile_pt.readline().strip()
	read = [line]
	read.append(seq)

	if in_format == "fastq":
		read.append(infile_pt.readline().strip())
		read.append(infile_pt.readline().strip())

	motifs = set()
	for x in xrange(min_motif,max_motif+1):
		motifs = motifs.union(set([seq[i:i+x] for i in xrange(0,len(seq),x)]))

	good = True
	for motif in motifs:
		new_seq = seq.replace(motif,'')
		if len(new_seq) <= (1 - threshold)*len(seq):
			good = False
			break

	if good:
		writeRead(read,outfile_good_pt,in_format)
		non_dusty += 1
	else:
		writeRead(read,outfile_bad_pt,in_format)
		dusty += 1
sys.stderr.write("Found " + str(non_dusty) + " non dusty reads and " + str(dusty) + " dusty reads\n")
