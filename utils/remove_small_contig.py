#!/usr/bin/env python3

from Bio import SeqIO
import sys

file=open(str(sys.argv[1])+"_sup200","w")

for sequence in SeqIO.parse(sys.argv[1], "fasta"):
	if len(sequence.seq) > 200 :
		file.write(">"+str(sequence.name) + "\n" + str(sequence.seq) + "\n")
