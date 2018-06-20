#!/usr/bin/python
import csv
import sys

csvin = csv.reader(sys.stdin,delimiter="\t")

row1 = csvin.next()
seq = []
while True:
	try:
		row = csvin.next()
		if row1[0] in row[0]:
			#row[1] = int(row1[1]) + int(row[1])
			row = [row[0]] + [int(row1[i]) + int(row[i]) for i in xrange(1,len(row))]
		else:
			seq.append(row1)
		row1 = row
	except StopIteration:
		seq.append(row1)
		break

csvout = csv.writer(sys.stdout,delimiter="\t")
csvout.writerows(seq)