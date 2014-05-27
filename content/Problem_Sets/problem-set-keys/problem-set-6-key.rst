.. _problem-set-6-key:

***********************
  Problem Set 6 : KEY
***********************

Problem 2.1
===========

The following R block calculates the range of each transcription factor
peak widths.

.. code-block:: r

    > library(dplyr)
    > colnames <- c('chrom','start','end','name')
    > peaks_df <- read.table('/opt//bio-workshop//data//peaks.bed.gz',
                              col.names=colnames)

    > peaks_tbl <- tbl_df(peaks_df)
    
    > summary <- peaks_tbl %.%
    +    group_by(name) %.%
    +    mutate(peak.width = end - start) %.%
    +    summarize(count = n(),
    +        min.width = min(peak.width),
    +        max.width = max(peak.width)) %.%
    +    mutate(range = max.width - min.width) %.%
    +    arrange(desc(range))

Now we want to grab the top ten factors ordered by range and make box and
violin plots.

.. code-block:: r

    > library(ggplot2)
    > factors.top10 <- summary$name[1:10]

    # the subset function returns parts of a data.frame, in this case it
    # returns rows where name is in the factors data.frame
    > peaks.subset <- subset(peaks_tbl, name %in% factors.top10)

    # end - start = width
    > gp <- ggplot(peaks.subset, aes(name, end - start, color=name))
    > gp + geom_boxplot() + scale_y_log10() + coord_flip()

    # save the plot
    ggsave('problem-2-1-boxplot.pdf')

    # or use geom_violin()
    > gp + geom_violin() + coord_flip()
    
    # save the plot
    ggsave('problem-2-1-violin.pdf')

Problem 2.2
===========

Grab the peaks that do not overlap any other peaks (besides itself). Strategy:

 #. Use :ref:`bedtools intersect <bedtools:intersect>`
 #. Give the same filename to ``-a`` and ``-b``
 #. Count the number of overlaps with ``-c``
 #. Use ``awk`` to select those regions with only 1 overlap (i.e. 1 means
    the only overlap is with itself)

.. code-block:: bash

    $ fname="peaks.bed.gz"
    $ bedtools intersect -a $fname -b $fname -c \
        | awk '$NF == 1' \
        | gzip -c \
        > peaks.no-overlap.bed.gz

Feed the new BED file into the pipeline for Problem 2.1 (above).

.. raw:: pdf

    PageBreak

