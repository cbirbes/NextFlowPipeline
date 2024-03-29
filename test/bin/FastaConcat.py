#!/bin/env python

import argparse
from Bio import SeqIO
import os

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--length', help = "chunck length", required=True)
args = vars(parser.parse_args())

chunck = args['length']

with open("pilonOut.fa","a") as output:
    FileIn = open("ctg.txt", "r")
    ListContigs=FileIn.read().split("\n")
    for x in range (0,len(ListContigs)):
        ListContigs[x]=ListContigs[x].replace(str(x+1)+";","").replace("\"","")
        for record in SeqIO.parse("pilonSR"+str(ListContigs[x])+".fasta","fasta"):
            if ("-" in ListContigs[x]):
                if (":1-" in ListContigs[x]):
                    output.write("\n>"+record.id+"\n")
                    output.write(str(record.seq))
                else:
                    output.write(str(record.seq))
            else:
                output.write("\n>"+record.id+"\n")
                output.write(str(record.seq))

cmd = "sed -i '1d' pilonOut.fa"
os.system(cmd)
