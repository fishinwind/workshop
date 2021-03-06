#! /usr/bin/env python

import sys
import pdb
from collections import defaultdict
from operator import itemgetter

if len(sys.argv) != 2:
    print >>sys.stderr, "error: specify bed file"
    sys.exit(1)

def parse_bed(bedfilename):
    ''' parse records and return each record '''

    for line in open(bedfilename):

        if line.startswith('#'):
            continue

        chrom, start, end, value = line.strip().split('\t')
        start = int(start)
        end = int(end)
        value = float(value)

        result = {'chrom':chrom, 'start':start, 'end':end, 'value':value}
        yield result

# for 1.
prev_start = None
# for 2.
coverage = defaultdict(int)

bedfilename = sys.argv[1]

for record in parse_bed(bedfilename):

    cur_start = record['start']

    # 1. calculate the distance between the start of the current record
    # and the previous record. note you may need to define a variable
    # outside of this loop.

    # this runs the *second* time through the loop, as prev_start is not
    # defined until a start is seen
    if prev_start:
        start_dist = cur_start - prev_start

        # XXX uncomment the next line to see the output

        # print "cur start dist is: %d" % start_dist

    # reset for next time through the loop
    prev_start = cur_start

    # 2. calculate the bases covered by the intervals for each
    # chromosome. note you may have to define a structure outside of
    # this loop to keep track of that information.

    cur_cov = record['end'] - record['start']
    coverage[record['chrom']] += cur_cov

# print out the coverage for each chromosome
for chrom in coverage:
    cov = coverage[chrom]
    print "coverage on %s = %s" % (chrom, cov)

# ------------------------------------------------------------------------

struct = defaultdict(list)

for record in parse_bed(bedfilename):
   
    chrom = record['chrom']
    
    # write additional code to get the start and end coordinates from
    # the record
    start = record['start']
    end = record['end']

    # create a tuple of coords. note where we put start and end:
    # coords[0] == start
    # coords[1] == end
    coords = (start, end)

    # add the coords to the growing list. replace `whichmeth` with the
    # appropriate method call

    # XXX whichmeth == append
    struct[chrom].append(coords)

for chrom in struct:
    # 1. use max() and min() in this loop to determine biggest start
    # values.

    # max() and min() look at the list of tuples, implicitly sort them by
    # size based on the first position, and then grab the min / max of the
    # sorted values. The `[0]` grabs the first item i.e. the start
    max_start = max(struct[chrom])[0]
    min_start = min(struct[chrom])[0]

    # 2. how do you change the max() and min() calls to look at the `end`
    # value instead of the `start`?

    # itemgetter is a function that you pass to the key argument to
    # fetch index 1 from the list, in this case index 1 = the `end` value.
    # the `[1]` grabs the second item i.e. the end

    max_end = max(struct[chrom], key=itemgetter(1))[1]
    min_end = min(struct[chrom], key=itemgetter(1))[1]

    # XXX: a different approach would be to flip the end and start
    # values above (line 80), and use max() / min() to look at those end
    # values

    # XXX: yet another approach, without itemgetter 
    # max_end = max(tup[1] for tup in struct[chrom])
    # min_end = min(tup[1] for tup in struct[chrom])

    print 'on chrom %s:' % chrom
    print '\tstarts: min = %s, max = %s' % (min_start, max_start)
    print '\tends: min = %s, max = %s' % (min_end, max_end)
