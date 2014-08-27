************
  Syllabus
************

General Information
===================

:Title:         Genome Informatics Workshop
:Course Number: MOLB FIXME 
:Semester:      Spring 2015
:Homepage:      http://hesselberthlab.github.io/workshop 
:Author:        Jay R. Hesselberth
:Author:        Brent Pedersen
:Organization:  University of Colorado School of Medicine
:Address:       Department of Biochemistry and Molecular Genetics
:Address:       Department of Medicine
:Copyright:     2013-2014 Jay R. Hesselberth, Brent Pedersen
:Copyright:     All Rights Reserved.
:Last updated:  |today|

PDF Content
===========

The course content is available as a combined PDF: 
:download:`PDF Download <_build/pdf/Genome-Informatics-Workshop.pdf>`

Instructor Information
======================

.. list-table::
    :widths: 40 40
    :header-rows: 1

    * - Instructor
      - Email
    * - Jay Hesselberth             
      - jay.hesselberth@gmail.com
    * - Brent Pedersen
      - bpederse@gmail.com
    * - XXX
      - XXX
    * - XXX
      - XXX

Schedule
========
Tue & Thurs, 1:00 - 3:00 PM. TA office hours from 3:30 PM
:ref:`See specific dates. <syllabus-specific-dates>`

Location
========
Health Sciences Library, Computer Teaching Lab 2

Course Description
==================

**The Genome Informatics Workshop is a hands-on tutorial of skills needed to
process large genomics data sets and visualize their results. The class
is taught from the standpoint of biologist with practical goals
(e.g. to interpret the results of a sequencing-based experiment and gain
biologically meaningful insight).**

We focus on working in the Linux environment, with emphasis on Linux
command-line tools, Python programming and the R statistical computing
environment. We use publicly available next-generation DNA sequencing data
from the ENCODE project to illustrate standard approaches for manipulating
sequencing data, aligning sequences to a reference genome, generating
coverage plots and displaying them in the UCSC Genome Browser. We will
cover specific analyses used in ENCODE project including ChIP-seq, DNase I
footprinting, mRNA-seq and genome sequencing to identify single nucleotide
and structural variants.

Course Credits
--------------

This is a 3 credit course.

Texts and Reading Materials
---------------------------

#. **Required**: Learn Python the Hard Way,
    http://learnpythonthehardway.org/book/

#. **Required**: Command Line Crash Course
    http://cli.learncodethehardway.org/book/

#. **Required**: ggplot2: Elegant Graphics for Data Analysis
    http://ggplot2.org/book/

#. **Required**: A Quick Guide to Organizing Computational Biology Projects
    http://www.ploscompbiol.org/article/metrics/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000424

Course Objectives
-----------------

  - Learn to manipulate large sequencing data sets with Linux tools
    and Python programming.

  - Learn to manipulate and visualize data in useful ways with the
    R statistical computing environment.

  - Learn workflows for experiments in ENCODE including ChIP-seq, DNaseI
    footprinting, mRNA-seq and variant detection.

  - Learn to visualize data from the ENCODE project in the UCSC Genome
    Browser

Canvas 
======

The course has a Canvas page [#]_ where announcements are made and
problem sets are uploaded

.. [#] https://ucdenver.instructure.com/courses/11079/assignments/syllabus

.. note::

    Be sure to **login** to the Canvas site to be able to see Announcements,
    upload Problem Sets, etc.

Assessment
==========

Progress of individual students will be assessed during the daily exercise
session, weekly problem sets, as well as a final project.

Grading Criteria
----------------

 - 50% participation
 - 40% problem sets (10 sets, 4% each)
 - 10% final project

.. _syllabus-specific-dates:

Specific Dates / Material to be Covered
=======================================

.. list-table::
    :widths: 20 40 80 20
    :header-rows: 1

    * - Class number
      - Date
      - Topic
      - Problem Set
    * - Class 1
      - T Jan 20
      - Introduction to VM, Linux and the shell
      - 
    * - Class 2 
      - Th Jan 22
      - Linux / Utilities
      - **PS1 due** (Mon Jan 26 12:00 AM)
    * - Class 3 
      - T Jan 27
      - Linux / Utilities
      - 
    * - Class 4 
      - Th Jan 29
      - Linux / Utilities
      - **PS2 due** (Mon Feb 212:00 AM)
    * - Class 5 
      - T Feb 3
      - Cluster Usage / Review
      - 
    * - Class 6 
      - Th Feb 5
      - Cluster Usage / Review
      - **PS3 due** (Mon Feb 9 12:00 AM)
    * - Class 7 
      - T Feb 10
      - Python
      - 
    * - Class 8 
      - Th Feb 12
      - Python
      - **PS4 due** (Mon Feb 16 12:00 AM)
    * - Class 9 
      - T Feb 17
      - Python 
      - 
    * - Class 10 
      - Th Feb 19
      - Python 
      - **PS5 due** (Mon Feb 23 12:00 AM)
    * - Class 11 
      - T Feb 24
      - Python 
      - 
    * - Class 12
      - Th Feb 26
      - ENCODE Overview
      - **PS6 due** (Mon Mar 2 12:00 AM)
    * - Class 13 
      - T Mar 3
      - BEDtools  
      - 
    * - Class 14 
      - Th Mar 5
      - ChIP-seq (coverage / peaks / motifs)
      - **PS6 due** (Mon Mar 9 12:00 AM)
    * - Class 15 
      - T Mar 10 
      - ChIP-seq / DNaseI-seq (UCSC)
      - 
    * - Class 16
      - Th Mar 12
      - Genomic analysis vignettes 
      - 
    * -
      - ** No Class Mar 16-20 (Campus Spring Break) **
      -
      - 
    * - Class 17 
      - T Mar 24
      - Genomic analysis vignettes 
      - 
    * - Class 18
      - Th Mar 26
      - R data & plotting 
      - **PS7 due** (Mon Mar 30 12:00 AM)
    * - Class 19
      - T Mar 31
      - R data & plotting 
      - 
    * - Class 20
      - Th Apr 2
      - R data & plotting 
      - **PS8 due** (Mon Apr 6 12:00 AM)
    * - Class 21
      - T Apr 7
      - R data & plotting 
      - 
    * - Class 22
      - Th Apr 9
      - R data & plotting 
      - **PS9 due** (Mon Apr 13 12:00 AM)
    * - Class 23
      - T Apr 14 
      - mRNA-seq (FPKM / diff exp)
      - 
    * - Class 24
      - Th Apr 16 
      - mRNA-seq (FPKM / diff exp)
      - **PS10 due** (Mon Apr 20 12:00 AM)
    * - Class 25 
      - T Apr 21
      - Exome Alignment
      - 
    * - Class 26 
      - Th Apr 23
      - Exome Variant Calling
      - 
    * - Class 27 
      - T Apr 28
      - TBD
      - 
    * - Class 28 
      - Th Apr 30
      - Final project presentations
      - 

.. raw:: pdf

    PageBreak
