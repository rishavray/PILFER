#!/usr/bin/python
import sys
import csv
import re

csvin = csv.reader(sys.stdin,delimiter="\t")
csvout = csv.writer(sys.stdout,delimiter="\t")
pattern = re.compile('NM:i:\d*')

for row in csvin:
	nm_id = [i for i, x in enumerate(row) if re.search(pattern, x)]
	nm_fields = row[nm_id[0]].split(':')
	ref_fields = row[2].split('|')

	#print str(len(ref_fields))
	edit_dist = int(nm_fields[2])
	seq_length = len(row[9])
	ref_seq_length = int(ref_fields[3])

	if (ref_seq_length - seq_length) <=2 and edit_dist <=1 :
		csvout.writerow(row)
