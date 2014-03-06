#! /usr/bin/env python

import sys
import pdb
import random

from itertools import combinations

from pybedtools import BedTool

if not len(sys.argv) >= 3:
    print >>sys.stderr, "error: specify at least 2 files to compare"
    sys.exit(1)

def get_cell_type(dirname):
    ''' take a directory, give me cell type '''
    fname = dirname.split('/')[-1].replace('wgEncodeUwDnase','')
    fname = fname.replace('PkRep1.narrowPeak.gz','')

    return fname

filenames = sys.argv[1:]
subset = random.sample(filenames, 10)

header_fields = ['#cell.type.1','cell.type.2','jaccard']
print '\t'.join(header_fields)

for fname1, fname2 in combinations(subset, r=2):

    tool1 = BedTool(fname1)
    tool2 = BedTool(fname2)

    result = tool1.jaccard(tool2, f=0.5, r=False)    
    stat = result['jaccard']

    ctype1 = get_cell_type(fname1)
    ctype2 = get_cell_type(fname2)

    fields = [ctype1, ctype2, stat]

    # print '\t'.join(fields)
    print '\t'.join(map(str, fields))

