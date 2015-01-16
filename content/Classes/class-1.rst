*****************************************************
             Class 1 : Class Introduction
*****************************************************

:Class date: Tues 2015 Jan 27

Goals
=====
#. Class overview
#. Get the VM running
#. Overall goals for the Blocks and Class

Class Overview
==============
Each class is 2 hours. We intend to spend the first 60 min going
through exercises that demonstrate how specific tools are useful in
bioinformatics. During the remaining hour, we expect you to work through
exercises, asking for help when you get stuck. 

We will record the first 60 minutes using Panopto Screen Capture, and
these recordings will be available in Canvas. We have found that simply
watching someone work in a terminal (move around, open up text editors,
write and execute simple programs) can be a very effective way to get
started with programming. Hopefully these movies will serve the same
purpose.

Each week, we will have 1 take home quiz, due the following Tuesday at 5
PM. 

Linux installations
===================
The PCs in the library have Virtual Box installed with a minimal Linux
installation. If you have your own PC laptop, you can install Virtual Box
and any standard Linux distribution (Ubuntu or Mint). If you have a Mac
laptop, you can do the same, or just use the native terminal.

In any case, we will create logins on our compute cluster (amc-tesla) and
all your work will be done through that.

Important Directories
=====================
You will learn about directories in the next class, but note that the
content for this class is available in:

.. code-block:: bash

    /opt/bio-workshop/

When you write answers to your homeworks, they should be runnable from
a directory like:

.. code-block:: bash

   /opt/bio-workshop/homeworks/section_1/week1/

Terminal and Gedit
==================
When you open a Terminal, you also launch a shell process, typically a
bash process. At the prompt, you can type things that bash understands,
and it will do them. The shell has its own language, which you will learn
over the course of the class. It also runs executable files that it can
find on the PATH (i.e. the set of directories that contain exectuables).

You can find what is on your PATH by typing:

.. code-block:: bash

   $ echo $PATH

The PATH is one of several environment variables that are created when you
login. You can see all of these with:

.. code-block:: bash

   $ env

.. nextslide::
    :increment:

One important program in the PATH is gedit. You will use this program to
keep notes and write small programs. There are other editors (e.g. vim and
emacs) that have additional features; it's up to you whether you would
like to learn those or stick with the graphical gedit.

To run gedit, use:

.. code-block:: bash

   $ gedit

This will open a window where you can type. You can write a small test
document and save it.

Shell and Python Programming
============================
It is important that you learn a few new computer languages. Others have
developed very good guides to teach you these languages, and we are going
to use those in the class. We expect you to begin taking these classes
immediately.

You will spend a lot of time going through these online classes, both in
scheduled class time, and outside of class time. Instead of focusing on
teaching you these languages, we will focus on helping you get through all
of the frustating problems that come up when you're learning the languages.

We will spend the first ~2 weeks learning shell [#]_ and all the things
you have access to within the shell.

.. [#] The Command Line Crash Course
        http://cli.learncodethehardway.org/book/

After learning the shell, we will begin learning R and several packages
within R.

Finally, we will begin learning Python [#]_. The Python language allows
you to do more sophisticated things that would be possible in shell or R, but
would be considerably more clunky.

.. [#] Learn Python the Hard Way
        http://learnpythonthehardway.org/book/

First Quiz : Reading
====================
Computational biology projects inevitably accrue a lot of files. For the
first quiz, you'll need to read a paper [#]_ and be able to put a set of
files in the correct places. We highly recommend adopting this scheme for
all of your projects in and out of the class.

.. [#] A Quick Guide to Organizing Computational Biology Projects (2009)
        PLoS Comput. Biol. William S. Noble
        http://dx.plos.org/10.1371/journal.pcbi.1000424

