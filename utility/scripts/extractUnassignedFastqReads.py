from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description='Extract unassigned fastq reads from a fastq file and an assignment file')
parser.add_argument('fq', type=str,help='input fastq')
parser.add_argument('a', type=str,help='input assignment file')
parser.add_argument('o', type=str,help='output fastq')

args = parser.parse_args()


with open(args.a,'r') as assignmentFile:
	lines = assignmentFile.read().splitlines()
	readIDs = []
	for line in lines:
		split = line.split('\t')
		readID = split[0]
		assignment = split[1]
		if assignment == '0':
			print("found unassigned read:{}".format(readID))
			readIDs.append(readID)
	reads = []
	for record in SeqIO.parse(args.fq,'fastq'):
		if record.id in readIDs:
			reads.append(record)
	SeqIO.write(reads,args.o,'fasta')
