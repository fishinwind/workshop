**************************
Class 2 : The command-line
**************************

:Class date: Thurs 29 January 2015

Goals
=====
#. The ``bash`` shell
#. continue learning to navigate within the terminal
#. understand the linux philosophy (small tools that do one thing well)
#. understand how to apply some common linux utilities to files
#. ``gedit`` to edit files

Unix Philosophy
===============
The Unix philosophy [#]_ Works well for bioinformatics:

+ Make each program do one thing well.
+ Make every program a filter.
+ Choose portability over efficiency.
+ Store data in flat text files.
+ Small is beautiful.
+ Use software leverage to your advantage.
+ Use shell scripts to increase leverage and portability.

.. [#] http://en.wikipedia.org/wiki/Unix_philosophy

Navigating In the Terminal
==========================
When you start the terminal, you will be in your home directory.

In Linux home is represented as `~` and also as `$HOME`

We will often show commands preceded with a '$' as you see in your
terminal.

Try this in the terminal:

.. code-block:: bash

    $ pwd

pwd stands for "print working directory"

.. nextslide::
    :increment:

Change to another directory:

.. code-block:: bash

    $ cd /tmp/

.. nextslide::
    :increment:

See what's in that directory:

.. code-block:: bash

    $ ls

Show more information:

.. code-block:: bash

    $ ls -lh

The `-lh` part are flags for the `ls` command.

These can also be separated like:

.. code-block:: bash

    $ ls -l -h

.. nextslide::
    :increment:

Getting Help In The Terminal
============================
How can you find out the arguments that `ls` accepts (or expects):

.. code-block:: bash

    $ man ls

and use spacebar to go through the pages. `man` is short for "manual" and
can be used on most of the commands that we will learn. 

In other linux software, it is common to get help by using:

.. code-block:: bash

    $ <program> -h

or:

.. code-block:: bash

    $ <program> --help

Which of these works for `ls`?

.. nextslide::
    :increment:

.. note::

    If you see an error message, read it carefully.  It may seem cryptic,
    but it is designed to inform you what went wrong.

Getting Help : Exercises
========================
+ use `man` to find out how to list files so that the most
  recently modified files are listed last.

(This is common when you're working on something and only
care about the most recently modified files)

+ use google to find the same thing. how else can you
  sort the output of `ls`?

Getting Help Outside The Terminal
=================================
*Use google*. Useful sites include:

+ stackexchange.com
+ biostars.org
+ seqanswers.com

In many cases, if you receive and error, you can copy-paste it into google
and find some info.

Other Commands In The Terminal
==============================
Use the `man` command to determine what `head` does.

Use `head` on the file /opt/bio-workshop/data/lamina.bed

Use `tail` to see the end of the file.

Questions
^^^^^^^^^
+ By default, head and tail show 10 lines. How can you see 13 lines?

Other Commands In The Terminal (Answers)
========================================

.. code-block:: bash

    $ man head
    $ head /opt/bio-workshop/data/lamina.bed
    $ tail /opt/bio-workshop/data/lamina.bed
    $ head -n 13 /opt/bio-workshop/data/lamina.bed

Word Counts (wc)
================
Exercise:

+ use **wc** to determine how many **lines** are in /opt/bio-workshop/data/lamina.bed
+ use **wc** to determine how many **words** are in /opt/bio-workshop/data/lamina.bed

Less (is More)
==============
To view a large file, use less::

    less /opt/bio-workshop/data/lamina.bed

You can forward-search in the file using "/"

You can backward-search in the file using "?"

You can see info about the file (including number of lines) using "ctrl+g"

You can exit **less** using "q"

Terminal History
================
Press the up arrow in the terminal.

Up and down arrows will allow you to scroll through your previous commands.

This is useful when running similar commands or when remembering what you have
done previously.

You can type the start of a command and then up-arrow and it will cycle
through commands that start with that prefix.

Tab-Completion
==============
The shell (bash) when set up properly can give you a lot of help

Type the following where [TAB] means the Tab key on the keyboard:

.. code-block:: bash

    $ cd /opt/bio-w[TAB]

Then hit tab. And:

.. code-block:: bash

    $ ls /opt/bio-w[TAB]

This will work for any file path and for any programs:

.. code-block:: bash

    $ hea[TAB]
 
.. nextslide::
    :increment:

What happens if you do:

.. code-block:: bash

    $ he[TAB][TAB] 

or:

.. code-block:: bash

    $ heaaa[TAB][TAB] 

Directory Shortcuts
===================
We have already used the `cd` command to change directories. And we have
used the "~" shortcut for home.

.. code-block:: bash

    $ cd ~ 
    $ ls ~

We can also move to or see what's in the parent directory with:

.. code-block:: bash
    
    $ ls ..
    $ cd ..

Or 3 directories up with:
    
.. code-block:: bash

    $ ls ../../..
    $ cd ../../..

To explicitly see the current directory:

.. code-block:: bash

    $ ls ./

.. nextslide::
    :increment:

We can go 2 directories up with:

.. code-block:: bash

    $ cd ../../

Here, we can remember that "." is the current directory and .. is one directory up.
What does this do:

.. code-block:: bash

    $ ls ./*

.. nextslide::
    :increment:

you can go to the last directory with:

.. code-block:: bash

    $ cd -

and switch back and forth by using that repeatedly.

Make and remove directories
===========================
.. code-block:: bash

    $ mkdir ~/tmp # OK

    $ mkdir ~/tmp/asdf/asdf # ERROR

    $ mkdir -p ~/tmp/asdf/asdf # OK


What does -p do?

Remove directories:

.. code-block:: bash

   $ rm ~/tmp/asdf # ERROR

   $ rm -r ~/tmp/asdf/asdf # OK

What does the -r flag do?

.. warning::

    Be careful with `rm -r` and `rm -rf`. You can accidentially remove
    entire directories that you didn't intend to.

Moving and copying files
========================
mv [source] [dest]:

.. code-block:: bash

    $ touch /tmp/asdf
    $ mv /tmp/asdf ~
    $ ls -lhtr ~/

.. nextslide::
    :increment:

In-class excercise:

#. make a directory `/tmp/moveable`
#. move that directory to ~
#. copy that directory to `/tmp/subdir/`

echo
====
`echo` means "print":

.. code-block:: bash

    $ echo "hello world"

and you can use it to see **bash** variables:

.. code-block:: bash

    $ echo $HOME

    $ echo $HISTFILE

Variables
=========
We will start covering programming in the next classes, but variables are a
key component of programming.

You can do::

    $ important=/opt/bio-workshop/data/lamina.bed
    $ ls -lh $important

sudo
====
.. image:: http://imgs.xkcd.com/comics/sandwich.png

.. code-block:: bash

    $ apt-get install cowsay
    $ sudo apt-get install cowsay

other commands
==============
excercise:

use `man` to determine the function of:

+ wget
+ uniq

How many records are present for each chromosome in
/opt/bio-workshop/data/lamina.bed (assume it is sorted by chromosome)?

gedit
=====
In order to edit files as you would using `notepad` or `word` in windows,
we will use the simple editor "gedit".

You can open gedit from the terminal using:

.. code-block:: bash

    $ gedit

This will open a new window with GUI controls. Use gedit to write/edit
scripts for this class.

Scripts
=======
A script is simply a series of commands that you save in a file. You will
need to write scripts to complete the homework.

Put this text:

.. code-block:: bash

    $ ls /opt/bio-workshop/

Into the file *`run.sh`* by opening `gedit` pasting that text then `save
as..` using the GUI controls

You can then run it as:

.. code-block:: bash

    $ bash run.sh

And you should see the same output as if you ran `ls /opt/bio-workshop` directly.

.. nextslide::
    :increment:

Scripts will be more useful when you have a series of commands you want to
run in series.

For example, a pipeline where you:

#. run quality control on some ChIP-seq reads 
#. align reads to a reference genome
#. find peaks (binding sites)
#. annotate the binding sites.

a script will provide a record of what you have done.

Scripts : Commenting
====================
For the homework, it if helpful to us if  you comment your scripts. 

Comments are not read by the shell, but they tell us (and you) what
you were trying to do. You can comment your code using the "#" symbol.

.. code-block:: bash
    
    # list all files in the /tmp/ directory ordered so that most recently
    # changed appear last
    $ ls -lhtr /tmp/

Pipes
=====
Since linux is made of small utilities, we often want to chain them
together. We will cover this in detail next class, but the idea
is that each program takes data, modifies it, and sends it to the next.

We can see lines 5-10 of a file with:

.. code-block:: bash

    $ head /opt/bio-workshop/data/lamina.bed | tail -n 5

.. raw:: pdf

    PageBreak
