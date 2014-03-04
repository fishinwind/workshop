#! /usr/bin/env python

import sys
import re

from itertools import combinations
from pyedtools import BedTool

filenames = sys.argv[1:]

def get_celltype(filename):
    regex = re.compile('')

if len(filenames) < 2:
    print >>sys.stderr, "error: specify >= 2 bed files"
    sys.exit(1)

for fname1, fname2 in combinations(filenames):
    
    btool1 = BedTool(fname1)
    btool2 = BedTool(fname2)

    stat = btool1.jaccard(btool2)

    cell_type_1 = get_celltype(fname1)
    cell_type_2 = get_celltype(fname2)

