
.. include:: /_static/substitutions.txt 

.. _problem-sets:

**************
 Problem Sets
**************

:Last updated: |today|

Problem Sets for the Genome Analysis Workshop (MOLB 7621):

.. toctree::
   :titlesonly:
   :glob: 

   Problem_Sets/problem_set_1/problem_set_1
   Problem_Sets/problem_set_2/problem_set_2
   Problem_Sets/problem_set_3/problem_set_3
   Problem_Sets/problem_set_4/problem_set_4
   Problem_Sets/problem_set_5/problem_set_5
   Problem_Sets/problem_set_6/problem_set_6
   Problem_Sets/problem_set_7/problem_set_7
   Problem_Sets/Problem_Set_Keys

Problem Set Keys
================

Past keys at :ref:`problem-set-keys`

.. _problem-set-submission:

Problem Set Submission
======================

In general, we want one run.sh file that includes all of the code
necessary to run the problem set.  This run.sh file should create any new
directories (ex. a dated directory in the results folder), perform
commands (ex. awk, cut, etc), and output results into a well-named file
(ex. > $results/$date/problem1.txt).  For each problem set, you should
also create a log file summarizing your results.

Once your problem set is complete, you will need to create a tar file.
Specify the root of your project directory and create a tar file of the
whole directory like this (change STUDENTID to your student ID):

.. code-block:: bash

    $ projectdir=$HOME/project
    $ tar -cvf STUDENTID-pset1.tar $projectdir

On the specific Problem Set assignment page at the Canvas site [#]_,
click the Submit Assignment button on the top right and paste the
full path to the tarfile (/vol3/home/username/.../foo.tar) in the text
box. Click the Submit Assignment button below the text box to complete 
the submission.

Auditors can e-mail their path to the TAs here: :ref:`contact-info`.

.. [#] https://ucdenver.instructure.com/courses/325063/assignments

.. raw:: pdf

    PageBreak

