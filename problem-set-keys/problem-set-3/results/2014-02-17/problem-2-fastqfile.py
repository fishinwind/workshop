#! /usr/bin/env python

'''
Problem 2.1 Which of the first 10 sequence records has the largest number of
'C' residues in the sequence? Report its record name.

Problem 2.2 For each of the first 10 records, Covert each character in the
quality score to a number, and sum the numbers. Use ord() to convert
characters to numbers

Problem 2.3
'''

import sys
from collections import Counter

if len(sys.argv) != 2:
    print >>sys.stderr, "error: specify fastq file"
    sys.exit(1)

fastqfilename = sys.argv[1]

def parse_fastq(fastqfile):
    ''' takes a file handle input and returns records as dicts '''

    name = seq = quals = None

    for line in fastqfile:

        if line.startswith('@'):
            # remove the @ char
            name = line.strip()[1:]
        elif line.strip() == '+':
            continue
        elif name and not quals and not seq:       
            seq = line.strip()
        else:
            quals = line.strip()
            yield {'name':name, 'seq':seq, 'quals':quals}
            name = seq = quals = None

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

print "Problem 2.1: the record name with the highest 'C' count is:"
print max(scores_2_1)[1]


print "Problem 2.2: sums of the first 10 quality scores in ord format:"
print '\n'.join(quals_2_2)

print "Problem 2.3:  most common sequences in the file:"
for seq, count in seqs.most_common(10):
    print seq,'\t',count

print "Problem 2.4: reverse complements of the first 10 seqs:"
print '\n'.join(quals_2_4)

