
.. include:: /_static/substitutions.txt

*************************
Class 18 : Applied Python 
*************************

:Class date: |c18-date|
:Last updated: |today|

Final Project
=============

We will have a final project that is worth 15% of your grade, due on April
2.

We want you to come up with a question you can address with the tools you
learned in the class. You should incorporate ENCODE data into your
analysis, and you can also integrate your own data.

We will go over some analysis vignettes on Wed and Fri to illustrate
example final projects.

Goals
=====

#. homework review
#. example application

Homework Review
===============

+ Everyone is having trouble with :ref:`problem-set-3`.

+ Trouble integrating simple python constructs (dicts, lists, strings) and
  logic (if blocks, for loops) to perform some tasks.

+ Specific questions about homework?

+ Questions 1, 2.

+ `toolshed module <https://pypi.python.org/pypi/toolshed>`_.

Application Impetus
===================

It is helpful to read python code that solves a common problem to
better understand how to use data-structures.

A python script that combines data from two files:

#. laboratory information on immune cell measurements.
#. information on a sequencing run for a subset of those samples.

This is real data (shuffled to protect innocent mice) and a real python script
(data from Ken Eyring and Ted Shade)

Data Files
==========

Download these files into "/vol1/opt/data/"

**laboratory info:** :download:`sample-lab-info.tsv </Miscellaneous/data/sample-lab-info.tsv>`

**sequence info:** :download:`sample-seq-info.csv </Miscellaneous/data/sample-seq-info.csv>`

**merging script:** :download:`sample-merge.py </Code/sample-merge.py>`

Once downloaded, look at the structure of the data files with `less`

We will spend this class *deriving*/*understanding* `sample-merge.py`

Set Up The Problem
==================

+ goal: merge the information from the 2 files.

+ we will add info from sample-seq-info.csv to sample-lab-info.tsv

+ sample-seq-info.csv contains a super-set of the samples in
  sample-lab-info.csv

+ we will match samples by the `Sample` column in sample-lab-info.tsv to
  the `Sample ID` column in sample-seq-info.csv

  * we will store rows from sample-seq-info.csv in a dictionary keyed by
    `Sample ID`

+ since we are **adding** to sample-lab-info.tsv, we don't need to filter
  out as we read from sample-seq-info.csv

Understand The Problem
======================

This is important!!

Any questions on file formats or our strategy?

Decide on Coding Strategy
=========================

#. read from `sample-seq-info.csv` into a dictionary with
   keys of `Sample ID` and values (dicts) with the data for that sample. We
   This dictionary is **seq_infos**. (skip header while fields[0] != "Lane")

#. loop (stream) over `sample-lab-info.tsv` to get **lab_info**
   for each sample

#. Find matching *sequence-info* for each row by using the `Sample` column as a
   key into **seq_infos**

#. The corresponding value of seq_infos[sample] will be all of the laboratory
   information for that sample.

#. Add the *seq_info* for the current sample to the *lab_info* using: 
   `lab_info.update(seq_info)`

#. print out the **lab_info** with newly added **seq_info**

Script
======

coming slides will go over the script block-by-block before viewing / 
running / modifying the entire script.

Script: Read seq info into dictionary
=====================================

.. code-block:: python

    # store data for all samples here, keys of sample-id, values of info
    seq_infos = {}

    # loop over each sample in seq_info
    for si in reader(seq_file, sep=",",
                     skip_while=is_extra_lines):
        sample_id = si['Sample ID']
        seq_infos[sample_id] = si

Now we have a dictionary with keys of sample ids and values of 
dictionaries containing the information for each sample.

We will use this as a lookup-table so that, given a sample_id from the
**lab_info** we can find the associated **seq_info**

Script: Iterate over lab-info and add seq-info
==============================================

We skip some error checking steps here for simplicity

.. code-block:: python

    is_first_line = True

    for lab_info in reader(lab_file):
        sample_id = lab_info['Sample']

        # we will add more logic here in the real script.
        seq_info = seq_infos[sample_id]
        lab_info.update(seq_info)
        # now lab_info has the sequene and the lab keys and values.

        if is_first_line: # print a header once only.
            print "\t".join(lab_info.keys()) 
            is_first_line = False

        # this will print out the data for each record.
        print "\t".join(lab_info.values()) 

Script: Run
===========

Let's run the script and see what comes out

.. code-block:: bash

    $ python example-merge.py > merged.tsv

look at merged info with `less` and verify that it has columns from
sample-lab-info.csv and sample-seq-info.csv


Script: Gedit
=============

Now let's open the script in gedit and go through it line-by-line!!

Script: Debug
=============

We can run the script from **ipython** as

.. ipython:: 

    In [1]: %run sample-merge.py

Open a gedit window and add some print statements to the script, followed by
"1/0" so that the script will stop and you can see what was printed. Save, then
run from ipython window.

This is a quick way to follow the flow of a script. As you understand each part,
move the print statement and the 1/0 further on in the script.
 
Spend the rest of class breaking, fixing and understanding this script.

.. raw:: pdf

    PageBreak
