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

    $ pwd

pwd is "print working directory"


Navigating In the Terminal (2)
==============================

Change to another directory:

    $ cd /tmp/

See what's in that directory:

    $ ls

Show more information:

    $ ls -lh

The "-lh" letters are flags, or modifier arguments to the *ls* command.
These can be separated as:

    $ ls -l -h

As you start working with a lot of files, you will want list them sorted
in order of modification so that the most recently modified appears last:

    $ ls -lhtr


Getting Help In The Terminal
============================

How can you find out the arguments that *ls* accepts (or expects)

    $ man ls

and use spacebar to go through the pages. *man* is short for manual
and can be used on all commands that we will learn. 

In other linux software, it is common to get help by using:

    $ program -h

or

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

    $ man head

    $ head ~/bio-workshop/data/some.fastq

    $ tail ~/bio-workshop/data/some.fastq

    $ head -n 13 ~/bio-workshop/data/some.fastq
        
    $ wc -l ~/bio-workshop/data/some.fastq



