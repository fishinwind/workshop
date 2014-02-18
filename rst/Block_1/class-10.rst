******************************
Class 10 : Intermediate Python 
******************************

Goals
=====

 #. homework review
 #. example application


Homework Review
===============

....

Application Impetus
===================

It is often helpful to read existing python code that solves a common
problem to better understand how to use certain data-structures.

We will go over a python script in class that combines data from two files:

#. laboratory information on immune cell measurements.
#. information on a sequencing run for a subset of those samples.


Data Files
==========

Download these files into "/opt/bio-workshop/data/"


**laboratory info:** :download:`sample-lab-info.tsv <../../data/sample-lab-info.tsv>`

**sequence info:** :download:`sample-seq-info.csv <../../data/sample-seq-info.csv>`

Once downloaded, look at the structure of these files with `less`


Set Up The Problem
==================

 + We will add info from sample-seq-info.csv to sample-lab-info.tsv

 + Note that sample-seq-info.csv contains a super-set of the samples in
   sample-lab-info.csv


 + we will match samples by the `Sample` column in sample-lab-info.tsv to
   the `Sample ID` column in sample-seq-info.csv

   * we will store rows from sample-seq-info.csv in a dictionary keyed by
     `Sample ID`

 + since we are adding to sample-lab-info.tsv, we don't even need to filter
   out as we read from sample-seq-info.csv

Decide on Coding Strategy
=========================

#. read info from `sample-lab-info.tsv` (tab-delimited) into a dictionary
   with keys of `Sample ID` and values of a dictionary containing all the
   data for that sample. We will call this dictionary **lab_info**

#. loop over `sample-seq-info` (comma-delimited) to get **sample_info** for
   each sample

   #. Find matching lab-info for each row by using the `Sample` column as a
      key into **lab_info**

   #. The corresponding value of lab_info[sample] will be all of the laboratory
      information for that sample.

   #. Add the lab_info for the current sample to the **sample_info** using:

    sample_info.update(lab_info)

   #. print out the **sample_info** with newly added **lab_info**


.. raw:: pdf

    PageBreak
