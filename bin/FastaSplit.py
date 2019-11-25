#!/bin/env python
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--length', help = "chuck length", required=True)
parser.add_argument('--fasta', help = "input fasta file", required=True)
args = vars(parser.parse_args())

FILENAME = args['fasta']
BPS_PER_FILE = args['length']  # base pairs per FASTA file
output = open("ctg.txt", "w")
compteur=1
y=0
contigsList=[]
totalLength=0

# There should be one and only one record, the entire genome:
for large_record in SeqIO.parse(FILENAME, "fasta"):
    length=len(large_record)
    totalLength+=int(length)

    if (int(length) < int(BPS_PER_FILE)):
        if (totalLength == int(length)):
            output.write(str(compteur)+";"+large_record.id)
        elif (totalLength < int(BPS_PER_FILE)):
            output.write(","+large_record.id)
        elif (totalLength == int(BPS_PER_FILE)):
            output.write(","+large_record.id+"\n")
        else:
            compteur+=1
            output.write("\n"+str(compteur)+";"+large_record.id)
            totalLength=length

    else:
        if (totalLength != length):
            output.write("\n")
            compteur+=1
        totalLength=0
        while (y+int(BPS_PER_FILE) <= length):
            output.write(str(compteur)+";"+large_record.id+':"'+str(y+1)+"-"+str(y+int(BPS_PER_FILE))+'"\n')
            compteur+=1
            y+=int(BPS_PER_FILE)
        if (y+int(BPS_PER_FILE) > length):
            output.write(str(compteur)+";"+large_record.id+':"'+str(y+1)+"-"+str(length)+'"\n')
            compteur+=1
    y=0

output.close()
