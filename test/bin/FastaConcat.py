#!/bin/env python

import argparse
import os
import os.path
from os import path
import glob

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--length', help = "chunck length", required=True)
args = vars(parser.parse_args())

chunck = args['length']

cmd = "sed -i '1d' *.fasta"
os.system(cmd)

y=1
x=1

items = os.listdir(".")
ListFile=[]
for names in items:
	if names.endswith(".fasta"):
		ListFile.append(names)
len1=len(ListFile)



while (x <= len1):
	y=1
	if (path.exists("pilonSRctg"+str(x)+".fasta")):
		y=1
		x+=1
	else:
		len1+=1
		with open("pilonSRctg"+str(x)+".fasta","a") as a1:
			while (glob.glob("pilonSRctg"+str(x)+":"+str(y)+"-*.fasta")):
				len1-=1
				with open(glob.glob("pilonSRctg"+str(x)+":"+str(y)+"-*.fasta")[0],"r") as a:
					for line in a:
						a1.write(line)
					a.close()
					y+=int(chunck)
			a1.close()
			x+=1

for x in range(1, len1+1):
	if (path.exists("pilonSRctg"+str(x)+".fasta")):
		with open("pilonSRctg"+str(x)+".fasta", 'r') as f:
			LenFile=0
			for line in f:
					LenFile+=len(line.replace("\n",""))
			f.close()
		with open("pilonSRctg"+str(x)+".fasta", 'r') as fx:
			with open("pilonOut.fa","a") as f1:
				if(LenFile > 0):
					f1.write("\n>ctg"+str(x)+" len="+str(LenFile)+"\n")
					for line in fx:
						f1.write(line.replace("\n",""))
			f1.close()
		fx.close()

cmd = "sed -i '1d' pilonOut.fa"
os.system(cmd)
