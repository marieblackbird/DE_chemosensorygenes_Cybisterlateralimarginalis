#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Python3 script to extract unique transcript IDs from a BLAST output

import os
import sys
import csv
from Bio import SeqIO

# Initialize a list to store unique transcript IDs
transcrit_odorant_receptor = []

# Open the BLAST output file provided as the first command-line argument
with open(sys.argv[1]) as csvfile:
    # Read the content of the file as a list of rows, each row being a list of values
    contenu_tableau = list(csv.reader(csvfile, delimiter='\t'))

# Loop through each row in the table
for line in contenu_tableau:
    # Check if the transcript ID (second column) is already in the list
    if line[1] not in transcrit_odorant_receptor:
        # If not, append it to the list
        transcrit_odorant_receptor.append(line[1])

# Print the number of unique transcripts found and the list of transcript IDs
print(str(len(transcrit_odorant_receptor)) + " transcript(s) candidate(s) \n" + str(transcrit_odorant_receptor))
