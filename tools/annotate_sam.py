#!/usr/bin/python
import csv
import sys

def ReverseComplement(seq):
	seq_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	return "".join([seq_dict[base] for base in reversed(seq)])

mini = 26
maxi = 33
seq_id = []
with open(sys.argv[1],"rb") as listin:
	for seq in listin:
		seq_id.append(seq.strip())

seq_id = list(set(seq_id))

csvin = csv.reader(sys.stdin, delimiter="\t")
csvout = csv.writer(sys.stdout, delimiter="\t")
for row in csvin:
	if not row[0][0] == "@":
		f = row[0].split(":")
		row[0] = f[1]
		seq = row[9]
		if (int(row[1]) & (0x10)):
			seq = ReverseComplement(seq)

		if seq in seq_id:
			row.append("XP:Z:PI")
		elif len(row[9])>=mini and len(row[9])<=maxi and int(f[0])>=100:
			row.append("XP:Z:PU")

		row.append("XC:i:"+f[0])
	csvout.writerow(row)

