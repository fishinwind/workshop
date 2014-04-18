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

cons_bedgraph = BedTool('/vol1/opt/data/hg19.100way.phyloP100way.bg.gz')
filenames = sys.argv[1:]

header_fields = ['#cell.type.1','cell.type.2','type','mean.cons']
print '\t'.join(header_fields)

for fname1, fname2 in combinations(filenames, r=2):

    tool1 = BedTool(fname1)
    tool2 = BedTool(fname2)

    # run this two times to get unique DHS from both comparisions
    result1 = tool1.intersect(tool2, v=True)    
    result2 = tool2.intersect(tool1, v=True)    

    cons1 = result1.map(cons_bedgraph, o='mean', c=4)
    cons2 = result2.map(cons_bedgraph, o='mean', c=4)

    pdb.set_trace()

    ctype1 = get_cell_type(fname1)
    ctype2 = get_cell_type(fname2)
    
    fields = [ctype1, ctype2, stat]

    # print '\t'.join(fields)
    print '\t'.join(map(str, fields))

