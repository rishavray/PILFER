#!/usr/bin/python
import csv
import sys 

#Use this script to merge a list of tsv files with a certain prefix, which are contained in a file passed as an argument
#The first argument is the file prefix, the second one is the base_dir

if len(sys.argv) < 3:
	sys.exit("Prefix list and/or base_dir not given")

prefix_list = []

with open(sys.argv[1],"rb") as listin:
	for prefix in listin:
		prefix_list.append(prefix.strip())

base_dir = sys.argv[2].rstrip("/")

union = {}
index = []

if len(sys.argv) == 4:
	arg = sys.argv[3]
	list_arg = arg.split(",")
	for ids in list_arg:
		if ':' in ids:
			temp = ids.split(':')
			index = index + [i for i in xrange(int(temp[0])-1,int(temp[1]))]
		else:
			index.append(int(ids)-1)
else:
	index = [i for i in xrange(0,len(prefix_list))]

#print index
i = 0
for idx in index:
	csvin = csv.reader(open(base_dir + "/"+prefix_list[idx],"rb"), delimiter="\t")
	
	for row in csvin:
		#print row[0]
		
		if row[0] in union:
			union[row[0]].append(row[1])
		else:
			union[row[0]] = [0 for x in xrange(0,i)]
			union[row[0]].append(row[1])

	for key in union:
		if len(union[key]) < i+1 :
			union[key].append(0)


	i += 1
#print union

csvout = csv.writer(sys.stdout,delimiter="\t")
csvout.writerow(["__"]+[prefix_list[i] for i in index])
for key in union:
	csvout.writerow([key]+union[key])

	#print key + "\t" + str(union[key])
