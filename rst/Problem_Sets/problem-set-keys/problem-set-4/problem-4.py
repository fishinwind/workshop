#! /usr/bin/env python

def parse_bed(bedfilename):
    ''' parse records and return each record '''

    for line in open(bedfilename):
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

bedfilename = 'lamina.bed'

for record in parse_bed(bedfilename):

    cur_start = record['start']

    # 1. calculate the distance between the start of the current record
    # and the previous record. note you may need to define a variable
    # outside of this loop.

    # this runs the *second* time through the loop, as prev_start is not
    # defined until a start is seen
    if prev_start:
        start_dist = cur_start - prev_start

    # XXX in your result, you would print this out. next line is commented
    # out so results file is not so cluttered

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

from collections import defaultdict
from operator import itemgetter

# specify the bedfilename 
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
    # sorted values
    max_start = max(struct[chrom])
    min_start = min(struct[chrom])

    # 2. how do you change the max() and min() calls to look at the `end`
    # value instead of the `start`?

    # itemgetter is a function that you pass to the key argument to
    # fetch index 1 from the list, in this case index 1 = the `end` value

    max_end = max(struct[chrom], key=itemgetter(1))    
    min_end = min(struct[chrom], key=itemgetter(1))    

    print 'on chrom %s:\n\t'
