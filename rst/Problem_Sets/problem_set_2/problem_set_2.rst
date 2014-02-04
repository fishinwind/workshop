Problem Set 2
=============

:Due date: 2014 Feb 18 at 5 PM MST

Overview
--------
For this quiz you will use several of the tools we learned about in the
last several classes, focusing on manipulating text files with Linux
command line tools.

.. note::

    Continue to use the organization scheme that we learned about in
    Problem Set 1. Part of our evaluation will include whether you are
    developing good organizational habits.

For all of the problems below. Your solution print out the region(s) as
they appear in the BED file with additional columns at the end. e.g.::

    chr12   1234    5678    0.93    my-extra-info

These should be separated by tabs (\\t), not spaces, unless otherwise indicated.

Remember you will combine `awk` with the other utilities you have learned.

Problem 1
---------

    What is the region with the largest start position (2nd column) on any
    chromosome in `lamina.bed`?

    What is the region with the largest end position on chrY in
    lamina.bed? Report this region in the format: "chr12:1234-5678"


Problem 2
---------

    What is the longest region (end - start) in `lamina.bed`?
    report as chrom\tstart\tend\tvalue\tregion_length

    What is the longest region (end - start) in `lamina.bed` with a value
    (4th column) greater than 0.9 on chr13. Report the header (1st line) in
    lamina.bed as well as the region.

Problem 3
---------

    What are the regions that overlap this interval in `lamina.bed`: 
    chr12:5,000,000-6,000,000? Report these so that they are ordered
    by descending value (4th column) and the columns are separated by commas
    rather than tabs.

    What is the average value (4th column) of those regions from `lamina.bed`
    that overlap the region: chr12:5,000,000-6,000,000?
    

Problem Set Submission
----------------------
Specify the root of your project directory and create a tar file of the whole
directory like this; you can change LASTNAME to your last name::

    $ projectdir=$HOME/bio-workshop/project
    $ tar -cvf LASTNAME-problem-set-2.tar $projectdir

Upload the tar file to the Problem Set at the Canvas site to complete the
submission.

