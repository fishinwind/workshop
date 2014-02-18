import sys
from toolshed import reader

seq_file = 'data/sample-seq-info.csv'
lab_file = 'data/sample-lab-info.tsv'

# this is a way to tell reader that we will skip all lines until we find one
# where fields[0] == "Lane"
def is_extra_lines(fields):
    return fields[0] != "Lane"

# we will store all the samples in seq_info, key'ed by the sample-id
seq_infos = {}
for si in reader(seq_file, sep=",",
                 skip_while=is_extra_lines):
    sample_id = si['Sample ID']
    seq_infos[sample_id] = si

# print this for 1 sample to see what it looks like:
print "!!! Key (sample id): %s\n!!! Value (dict): %s\n" % seq_infos.iteritems().next()

# now we can use that dictionary as a lookup while we iterate through the
# lab-info file.

not_found = []

header_printed = False
for line_no, lab_info in enumerate(reader(lab_file)):

    sample_id = lab_info['Sample']
    # now the the seq info associated with that sample
    try:
        seq_info = seq_infos[sample_id]
    except KeyError:
        not_found.append(sample_id)
        continue

    # now update the current lab_info dict with the seq_info
    lab_info.update(seq_info)

    # if it is the first line
    if not header_printed:
        print "\t".join(lab_info.keys())

    print "\t".join(lab_info.values())

print "samples not found:", ",".join(not_found)
