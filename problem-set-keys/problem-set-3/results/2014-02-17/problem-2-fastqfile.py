#! /usr/bin/env python

'''Which of the first 10 sequence records has the largest number of
'C' residues in the sequence? Report its record name'''

import sys
from collections import Counter
import pdb

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

seen_records = 0
record_limit = 10
scores = []

for record in parse_fastq(open(fastqfilename)):

    nuc_counts = Counter(record['seq'])

    # add a tuple of (C_count, name) to scores
    scores.append((nuc_counts['C'], record['name']))

    seen_records += 1

    if seen_records == record_limit:
        break

print "the record name with the highest 'C' count is:"
print max(scores)[1]
