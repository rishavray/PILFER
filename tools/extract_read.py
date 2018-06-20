#!/usr/bin/python
import sys
import csv

def ReverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])

csvin = csv.reader(sys.stdin,delimiter="\t")

for row in csvin:
	count = int(row[len(row)-1][5:])
	tag = row[len(row)-2][5:]
	if tag =="PI" or tag == "PU":
		seq = row[9]
		strand = "+"
		if (int(row[1]) & (0x10)):
			seq = ReverseComplement(seq)
			strand = "-"
		
		print row[2] + "\t" + row[3] + "\t" + str(int(row[3])+len(seq)-1) + "\t" + seq + "::" + tag + "\t" + str(count)  + "\t" + strand

