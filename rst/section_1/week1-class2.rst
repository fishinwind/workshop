Class 2: The command-line
=========================

Goals
=====

1. Introduce IPython notebook
2. continue learning to navigate within the terminal
3. understand the linux philosophy (small tools that do one thing well)
4. understand how to apply some common linux utilities to files
5. understand how to chain/stream/pipe linux operations


Unix Philosophy
===============

+ Small is beautiful.
+ Make each program do one thing well.
+ Choose portability over efficiency.
+ Store data in flat text files.
+ Use software leverage to your advantage.
+ Use shell scripts to increase leverage and portability.
+ Make every program a filter.

(From Mike Gancarz http://en.wikipedia.org/wiki/Unix_philosophy)


Navigating In the Terminal
==========================

When you start the terminal, you will be in your home directory.

In Linux home is represented as ~

We will show commands preceded with a '$' as you see in your terminal

Try this in the terminal:

.. code-block:: bash

    $ pwd

pwd is "print working directory"


Navigating In the Terminal (2)
==============================

Change to another directory:

.. code-block:: bash

    $ cd /tmp/

See what's in that directory:

.. code-block:: bash

    $ ls

Show more information:

.. code-block:: bash

    $ ls -lh

The "-lh" letters are flags, or modifier arguments to the *ls* command.
These can be separated as:

.. code-block:: bash

    $ ls -l -h

As you start working with a lot of files, you will want list them sorted
in order of modification so that the most recently modified appears last:

.. code-block:: bash

    $ ls -lhtr


Getting Help In The Terminal
============================

How can you find out the arguments that *ls* accepts (or expects)

.. code-block:: bash

    $ man ls

and use spacebar to go through the pages. *man* is short for manual
and can be used on all commands that we will learn. 

In other linux software, it is common to get help by using:

.. code-block:: bash

    $ program -h

or

.. code-block:: bash

    $ program --help



Getting Help Outside The Terminal
=================================

Use google. Favor results on:

 + stackexchange.com
 + biostars.org
 + seqanswers.com

In many cases, if you receive and error, you can copy-paste it into google and find some info.


Other Commands In The Terminal
==============================

Use the *man* command to determine what *head* does.

Use *head* on the file ~/bio-workshop/data/some.fastq

Use *tail* to see the end of the file.

By default, head and tail show 10 lines. How can you see 13 lines?

How many lines are in the file. Use *wc*


Other Commands In The Terminal (Answers)
========================================

.. code-block:: bash

    $ man head

    $ head ~/bio-workshop/data/some.fastq

    $ tail ~/bio-workshop/data/some.fastq

    $ head -n 13 ~/bio-workshop/data/some.fastq
        
    $ wc -l ~/bio-workshop/data/some.fastq


Terminal History
================

Press the up arrow in the terminal.

Up and down arrows will allow you to scroll through your previous commands.

This is useful when running similar commands or when remembering what you have
done previously.


Tab-Completion
==============

The shell (bash) when set up properly can give you a lot of help
Type the following where [TAB] means the Tab key on the keyboard:

.. code-block:: bash

    $ cd ~/bio-w[TAB]

Then hit tab. And:

.. code-block:: bash

    $ ls ~/bio-w[TAB]

This will work for any file path.


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

We can go 2 directories up with:

.. code-block:: bash

    $ cd ../../

Here, we can remember that "." is the current directory and .. is one directory up.
What does this do:

.. code-block:: bash

    $ ls ./*

Scripts
=======

A script is simply a series of commands that you save in a file. You will need to write
scripts to complete the homework.

Put this text:

    ls ~/bio-workshop/

Into the file `my-ls.sh`.

You can then run it as:

.. code-block:: bash

    bash my-ls.sh

And you should see the same output as if you ran `ls ~/bio-workshop` directly.

Scripts
=======

Scripts will be more useful when you have a series of commands you want to run in series.

For example a pipeline where you:

 1. run quality control on some ChIP-seq reads 
 2. align reads to a reference genome
 3. find peaks (binding sites)
 4. annotate the binding sites.

In cases like that, a script will provide a record of what you have done.

