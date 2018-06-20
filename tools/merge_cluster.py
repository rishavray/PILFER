#!/usr/bin/python
import csv
import sys

bed_file = open(sys.argv[1],"rb")
cluster_file = open(sys.argv[2],"rb")
size = int(sys.argv[3])
clusters = []
empty = [0 for x in xrange(0,size)]

for line in bed_file:
	line = line.strip()
	clusters.append(line.split() + empty)

csvin = csv.reader(cluster_file,delimiter="\t")
next(csvin)
for row in csvin:
	field = row[0].split(":")
	#print field
	pos = field[1].split("-")
	pos = [int(x) for x in pos]

	for clust in clusters:
		if field[0] == clust[0] and int(clust[1]) <= pos[0] and int(clust[2]) >= pos[1]:
			for i in xrange(3,size+3):
				clust[i] = float(clust[i]) + float(row[i-2])

csvout = csv.writer(sys.stdout,delimiter="\t")
for row in clusters:
	row[0] = row[0] + ":" + str(row[1]) + "-" + str(row[2])
	del row[1:3]
csvout.writerows(clusters)
