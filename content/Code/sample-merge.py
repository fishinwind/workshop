import sys
from toolshed import reader

# define our input files
seq_file = 'sample-seq-info.csv'
lab_file = 'sample-lab-info.tsv'


# this is a way to tell reader that we will skip all lines until we find one
# where fields[0] == "Lane"
def is_extra_lines(fields):
    return fields[0] != "Lane"

# we will store all the samples in seq_info, key'ed by the sample-id
seq_infos = {}

# loop over each sample in seq_info
for si in reader(seq_file, sep=",",
                 skip_while=is_extra_lines):
    sample_id = si['Sample ID']
    # NOTE: print si here to see what it looks like.
    seq_infos[sample_id] = si

# NOTE: print this for 1 sample to see what it looks like:
# print "!!! Key (sample id): %s\n!!! Value (dict): %s\n" % seq_infos.iteritems().next()

# NOTE: print the entire dictionary of dictionaries:
# print seq_infos

# NOTE: stop here until you understand what's going on above!!!

# now we can use that dictionary as a lookup while we iterate through the
# lab-info file.

# We will store the id's in the lab_info that we didn't find in the seq_info
# so we can tell the researcher about the missing seq data.
not_found = []

is_first_line = True
for lab_info in reader(lab_file):
    # NOTE: print lab_info here to see what it looks like

    sample_id = lab_info['Sample']
    # now the the seq info associated with that sample
    try:
        seq_info = seq_infos[sample_id]
    except KeyError:
        not_found.append(sample_id)
        continue

    # now update the current lab_info dict with the seq_info
    lab_info.update(seq_info)
    # NOTE: print lab_info here to see what it looks like now that it
    #       contains the seq_info as well.



    if is_first_line:
        # NOTE: order is not preserved in a dictionary so the header will have
        #       the seq columns interspersed with the lab columns. But, the
        #       keys and values will always come out in the same order between
        #       rows.
        print "\t".join(lab_info.keys())
        # then set to true so we only print the header once.
        is_first_line = False

    # print out the values for every record.
    print "\t".join(lab_info.values())

print >>sys.stderr, "samples not found in seq data:", ",".join(not_found)
