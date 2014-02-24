#! /usr/bin/env python

'''
Problem Set 3 (FASTQ)
=====================

Problem 2.1
-----------
Which of the first 10 sequence records has the largest number of 'C'
residues in the sequence? Report its record name.

Problem 2.2
-----------
For each of the first 10 records, Covert each character in the quality
score to a number, and sum the numbers. Use ord() to convert characters to
numbers.

Problem 2.3
-----------
Use the Python Counter to count unique sequences in the Fastq file. Report
the top ten most abundant sequences in the Fastq file

Problem 2.4
-----------
Report the revese complement of each of the first 10 sequences.
'''

import sys
from collections import Counter

if len(sys.argv) != 2:
    print >>sys.stderr, "error: specify fastq file"
    sys.exit(1)

fastqfilename = sys.argv[1]

# --- define functions -------------------------------------------
def parse_fastq(fastqfile):
    ''' takes a file handle input and returns records as dicts '''
    line_num = 0  # Keep track of lines in file
    
    for line in fastqfile:

        if line_num % 4 == 0:
            # Removes @ character from record name
            name = line.strip()[1:]  

        elif line_num % 4 == 1:      
            seq = line.strip()

        elif line_num % 4 == 3:
            quals = line.strip()
            yield {'name':name, 'seq':seq, 'quals':quals}

        # line incremented with every iteration of loop
        line_num += 1

def reverse_comp(seq):
    ''' reverse complement a sequence. return in 5'->3' direction '''
    rc_seq = []
    for char in reversed(seq):
        if char == 'T':
            rc_seq.append('A')
        elif char == 'A':
            rc_seq.append('T')
        elif char == 'G':
            rc_seq.append('C')
        if char == 'C':
            rc_seq.append('G')

    return ''.join(rc_seq)

# --- main program -------------------------------------------

# keep track of record number for first 10
seen_records = 0
record_limit = 10

# problem 2.1 data
scores_2_1 = []

# problem 2.2 data
quals_2_2 = []

# problem 2.4 data
quals_2_4 = []

for record in parse_fastq(open(fastqfilename)):

    nuc_counts = Counter(record['seq'])

    # add a tuple of (C_count, name) to scores
    scores_2_1.append((nuc_counts['C'], record['name']))

    # problem 2.2 data
    qual_sum = sum([ord(char) for char in record['quals']])
    # covert to str for printing later
    quals_2_2.append(str(qual_sum))
   
    quals_2_4.append(reverse_comp(record['seq']))

    seen_records += 1

    if seen_records == record_limit:
        break

# problem 2.3 data
seqs = Counter()
for record in parse_fastq(open(fastqfilename)):
    seqs[record['seq']] += 1 

# --- print results ---------------------------------------------------

print "Problem 2.1: the record name with the highest 'C' count is:\n"
print max(scores_2_1)[1]

# blank line
print

print "Problem 2.2: sums of the first 10 quality scores in ord format:\n"
print '\n'.join(quals_2_2)

# blank line
print

print "Problem 2.3:  most common sequences in the file:"
for seq, count in seqs.most_common(10):
    print seq,'\t',count

# blank line
print

print "Problem 2.4: reverse complements of the first 10 seqs:\n"
print '\n'.join(quals_2_4)

