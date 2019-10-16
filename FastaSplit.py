import argparse

parser = argparse.ArgumentParser(description='Description of your program')
parser.add_argument('--length', help = "chuck length", required=True)
parser.add_argument('--fasta', help = "input fasta file", required=True)
args = vars(parser.parse_args())

FILENAME = args['fasta']

BPS_PER_FILE = args['length']  # base pairs per FASTA file

# There should be one and only one record, the entire genome:

output = open("ctg.txt", "w")

print(FILENAME)
large_record = open(FILENAME, 'r')
large_record = large_record.read()
large_record = large_record.replace("A\n","").replace("C\n","").replace("G\n","").replace("T\n","")
large_record = large_record.replace("A>","A\n>").replace("C>","C\n>").replace("G>","G\n>").replace("T>","T\n>")
large_record = large_record.replace(" ","\n")
large_record = large_record.split("\n")
y = len(large_record)
x=0
compteur=1
while (x < y):
	if (large_record[x][0]!=">"):
		print("Erreur fichier fasta")
	if ("RC" in large_record[x+2]):
		del large_record[x+2]
	if ("XC" in large_record[x+2]):
		del large_record[x+2]
	nameContig = large_record[x].replace(">","")
	# ~ print (nameContig)
	lenContig = large_record[x+1].replace("len=","")
	lenContig = lenContig.replace("LN:i:","")
	# ~ print(lenContig)
	if (int(lenContig) < int(BPS_PER_FILE)):
		output.write(str(compteur)+";"+nameContig+"\n")
		compteur+=1
	else:
		for y in range (0, int(lenContig), int(BPS_PER_FILE)):
			if ((y + int(BPS_PER_FILE)) > int(lenContig)):
				output.write(str(compteur)+";"+nameContig+':"'+str(y+1)+"-"+lenContig+'"\n')
				compteur+=1
			else:
				z=y+int(BPS_PER_FILE)
				output.write(str(compteur)+";"+nameContig+':"'+str(y+1)+"-"+str(z)+'"\n')
				compteur+=1
	y = len(large_record)
	x+=3
output.close()
