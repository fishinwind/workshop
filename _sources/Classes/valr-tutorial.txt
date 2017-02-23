
.. include:: /_static/substitutions.txt

=============================
Class 11 : Genome Arthimetic with valr 
=============================

:Last updated: |today|

Goals
-----
#. Introduce the valr package for performing genome arithmetic in R
#. Review metagene analyses 
#. Introduce randomizing intervals

.. |valr| image:: ../_static/images/valr_logo.png

|valr| philosophy behind valr
---------------------

#. Promote fluid interactive analysis by limiting back-and-forth between
   command line and R. 
#. Empower data parsing/manipulating power of ``dplyr`` and
   the ``tidyverse`` style of programming instead of custom data formats
#. Encourage the use of RMarkdown and Shiny apps for reproducible and
   accessible genomics analysis  

Installation
------------

.. code-block:: r
    
    install.packages('valr')
    library('valr')

Documentation
-------------

The ``valr`` package has extensive documentation that can be accessed via
R's help interface

.. code-block:: r
   
    ?valr
    ?read_bed
    ?bed_intersect
    ?sorting
    
Alternatively, valr's `online documentation
<http://rnabioco.github.io/valr>`_ contains documentation for each
function as well as examples and vignettes.  

Basic Usage
-----------

#. read bed/bedgraph/vcf data into R as tbls
   (``read_bed()``, ``read_genome()``)
#. all bed functions start with ``bed_*()``
   (``bed_intersect()``, ``bed_map()``)
#. Pass data via ``%>%`` 
#. Execute arbitrary aggregating functions in summarizing operations
   (``bed_map(x, y, my_custom_function(score)``)

Example
-------

.. code-block:: r
    
    library(valr)
    library(dplyr)
    
    snps <- read_bed(valr_example('hg19.snps147.chr22.bed.gz'), n_fields = 6)
    genes <- read_bed(valr_example('genes.hg19.chr22.bed.gz'), n_fields = 6)
    
    # find snps in intergenic regions
    intergenic <- bed_subtract(snps, genes)
    # find distance from intergenic snps to nearest gene
    nearby <- bed_closest(intergenic, genes)
    
    nearby %>%
      select(starts_with('name'), .overlap, .dist) %>%
        filter(abs(.dist) < 5000)
    
Grouping 
-------

valr functions respect groupings set by ``dplyr::group_by()``

.. code-block:: r
   
    # to intersect by strand
    snps <- group_by(snps, strand)
    genes <- group_by(genes, strand)
    bed_intersect(snps, genes)

Summaries by Column
-------------------

``valr`` functions accept named columns and permit mutiple name/value
summaries 

.. code-block:: bash
   
    # calculate the mean of column 6 for intervals in `b` that overlap
    # with `a`
    bedtools map -a a.bed -b b.bed -c 6 -o mean

.. code-block:: r
    
    # calculate the mean and variance for a `value` column
    bed_map(a, b, .mean = mean(value), .var = var(value))
    
    # report concatenated and max values for merged intervals
    bed_merge(a, .concat = concat(value), .max = max(value))

Random functions
----------------

``valr`` and ``bedtools`` have a series of functions useful for generating
background interval datasets for statistical tests

.. code-block:: r
    
    bed_random(genome) # random intervals of a fixed sized
    bed_shuffle(bed, genome) # randomly placed intervals of same input size
    bed_flank(bed, both = 100) # get flanking regions
    dplyr::sample_n() # get random rows 

Fast Computation
----------------

Computationally intensive functions in ``valr`` are written in Rcpp/C++,
making ``valr`` fast enough for interactive analysis

.. code-block:: r

    genomefile <- valr_example('hg19.chrom.sizes.gz')
    x <- bed_random(genome) # generate 1e6 random 1kp ivls
    y <- bed_random(genome)
    bed_intersect(x, y) #should take ~1 second

TSS Metagene
------------

Previously we used bedtools to examine the coverage of
CTCF ChIP-Seq data over Transcription Start Sites. 

A similar analysis can easily be performed with ``dplyr`` and ``valr``

See a detailed example using valr on `H3K4me3 ChIP-Seq
<http://rnabioco.github.io/valr/articles/valr.html#meta-analysis>`_

Load data
---------

.. code-block:: r

    # `valr_example()` identifies the path of example files
    bedfile <- valr_example('genes.hg19.chr22.bed.gz')
    genomefile <- valr_example('hg19.chrom.sizes.gz')
    bgfile  <- '~/data-sets/bedtools/ctcf.hela.chr22.bg.gz'

    genes <- read_bed(bedfile, n_fields = 6)
    genome <- read_genome(genomefile)
    y <- read_bedgraph(bgfile)


Make TSS intervals
------------------

.. code-block:: r
    
    # generate 1 bp TSS intervals, `+` strand only
    tss <- genes %>%
      filter(strand == '+') %>%
      mutate(end = start + 1)

    region_size <- 1000 # 1000 bp up and downstream
    win_size <- 5  # 5 bp windows

    # add slop to the TSS, break into windows and add a group
    x <- tss %>%
      bed_slop(genome, both = region_size) %>%
      bed_makewindows(genome, win_size)
    x
    #> # A tibble: 13,530 × 7
    #>    chrom    start      end      name score strand .win_id
    #>    <chr>    <dbl>    <dbl>     <chr> <chr>  <chr>   <dbl>
    #> 1  chr22 16161065 16161115 LINC00516     3      +       1
    #> 2  chr22 16161115 16161165 LINC00516     3      +       2
    #> 3  chr22 16161165 16161215 LINC00516     3      +       3

Compute and Summarize Coverage
------------------------------

.. code-block:: r
    
    # map signals to TSS regions and calculate summary statistics.
    res <- bed_map(x, y, win_mean = mean(value, na.rm = TRUE)) %>%
      group_by(.win_id) %>%
      summarize(win_sum = sum(win_mean, na.rm = TRUE))

    res
    
    #> # A tibble: 401 × 2
    #>   .win_id win_sum
    #>     <dbl>   <dbl>
    #> 1        1   205.0
    #> 2        2   203.5
    #> 3        3   207.0
    #> 4        4   215.0
    #> 5        5   222.5

Plot metagene
-------------

.. code-block:: r
    
    library(ggplot2)
    
    x_labels <- seq(-region_size, region_size, by = win_size * 50)
    x_breaks <- seq(1, 401, by = 50)
    
    ggplot(res, aes(x = .win_id, y = win_sum)) +
      geom_point() + 
      scale_x_continuous(labels = x_labels, breaks = x_breaks) +
      xlab('Position (bp from TSS)') + ylab('Signal') +
      ggtitle('Human CTCF signal near transcription start sites') +
      theme_classic()


Exercises
--------

#. Use ``valr`` to generate a TSS metagene plot
#. Add a second metagene plot line from randomized or shuffled regions (hint: randomize TSS positions)
#. Plot metagene with TSS's derivied from both positive and negative
   strands (hint: reverse negative strand interval windows)

Contributing
------------

If you are interested please contribute feature requests, new code,
documentation or your ideas for analysis vignettes.

https://github.com/rnabioco/valr

