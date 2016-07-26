#!/usr/bin/env python2.7

#Created by: Bruno Costa	
# ITQB 2016
# 
# This calculates the distribution profile from collapsed fasta files
# fragment-abundance.py -h for help

import argparse

##### Inputs #######
#targets_file="/home/brunocosta/INIA-Targets/INIA-miRNA-classification_INIA_mappable.fa_dd__Cleaveland4_result.tsv"
#write_targets="Ppi-INIA-sRNA-all-Cleavelane4-Results-FragCount.tsv"
#density_file="/home/brunocosta/INIA-Targets/INIA_mappable.fa_dd.txt"
#

parser = argparse.ArgumentParser(description='This program add the fragment count to the cleaveland 4 output tsv based on the degradome density file')

parser.add_argument('--targets',type=str, nargs=1, metavar='Targets file',dest='targets',required=True,help='Path to the file outputed by cleaveland 4 with headers and metadata')
parser.add_argument('--output',type=str, nargs=1, metavar='Output file',dest='output',required=True,help='The output file name and destination. ex: /home/user/outputfile.tsv')
parser.add_argument('--density',type=str, nargs=1, metavar='Degradome Density file',dest='density',required=True,help='The path to the degradome density file')

args = parser.parse_args()

targets_file=args.targets[0]
write_targets=args.output[0]
density_file=args.density[0]





targets=open(targets_file,"r")
density=open(density_file,"r")

writer=open(write_targets,"w")
#This is a script provided by FPMartins to index fasta into a dictionary	
d= {}
#reading from density file
for lines in density:
	if lines.startswith('@ID'):
		name = lines.strip().split(":",1)[1]
		d[name]={}
	else:
		if lines.startswith("@LN"):
			length_of_unigene=lines[1:].split(":",1)
		else:
			#Append the line to the dictionary inside the dictionary for this entry
			line=lines.strip().split("\t")
			if (len(line))==3:
				d[name][line[0]]=line[1:] 

density.close()


#Parses the target file excluding the metadata information in the beginging.
annot_targets=[lines.strip().split("\t") for lines in targets if not lines.startswith("#") ]

header=annot_targets[0]

res=[]
#Skip header line
for annotation in annot_targets[1:]:
	#Grab fragment info about transcript and tslice site 
	annotation.append(d[annotation[2]][annotation[5]][0])
	res.append(annotation)

#Output to file
writer.write(reduce(lambda x,y:x+"\t"+y, header)+"\tFrag #\n")
for i in res:	
	output=reduce(lambda x,y:x+"\t"+y, i)+"\n"
	writer.write(output)
writer.flush()	
writer.close()
