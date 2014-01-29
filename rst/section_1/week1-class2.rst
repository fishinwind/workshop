Class 2: The command-line
=========================

Goals
-----

1. The **Shell**
2. continue learning to navigate within the terminal
3. understand the linux philosophy (small tools that do one thing well)
4. understand how to apply some common linux utilities to files
5. gedit to edit files


Unix Philosophy
---------------

Works well for bioinformatics:

+ Make each program do one thing well.
+ Make every program a filter.
+ Choose portability over efficiency.
+ Store data in flat text files.
+ Small is beautiful.
+ Use software leverage to your advantage.
+ Use shell scripts to increase leverage and portability.

(From Mike Gancarz http://en.wikipedia.org/wiki/Unix_philosophy)


Navigating In the Terminal
--------------------------

When you start the terminal, you will be in your home directory.

In Linux home is represented as `~` and also as `$HOME`

We will often show commands preceded with a '$' as you see in your terminal

Try this in the terminal:

.. code-block:: bash

    $ pwd

pwd is "print working directory"


Navigating In the Terminal (2)
------------------------------

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

Navigating In the Terminal (3)
------------------------------

The content for the class, including the source for generating this document,
are in:

.. code-block:: bash

    /opt/bio-workshop/

on your virtual machine.

Navigate to that directory and look around with `ls`


Getting Help In The Terminal
----------------------------

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

Which of those works for `ls`?

Getting Help: Exercises
-----------------------


 + use `man` to find out how to list files so that the most
   recently modified files are listed last.

(This is common when you're working on something and only
care about the most recently modified files)

 + use google to find the same thing. how else can you
   sort the output of `ls`?


Getting Help In The Terminal(2)
-------------------------------

 + If you see an error message, read it carefully. 
 + It may seem cryptic, but it is built to inform you what went wrong.


Getting Help Outside The Terminal
---------------------------------

Use google. Favor results on:

 + stackexchange.com
 + biostars.org
 + seqanswers.com

In many cases, if you receive and error, you can copy-paste it into google and find some info.


Other Commands In The Terminal
------------------------------

Use the *man* command to determine what *head* does.

Use *head* on the file /opt/bio-workshop/data/lamina.bed

Use *tail* to see the end of the file.

Questions
+++++++++

By default, head and tail show 10 lines. How can you see 13 lines?

How many lines are in the file. (Use *wc*)


Other Commands In The Terminal (Answers)
----------------------------------------

.. code-block:: bash

    $ man head

    $ head /opt/bio-workshop/data/lamina.bed

    $ tail /opt/bio-workshop/data/lamina.bed

    $ head -n 13 /opt/bio-workshop/data/lamina.bed
        
    $ wc -l /opt/bio-workshop/data/lamina.bed

Word Counts
-----------

Exercise:

    + use **wc** to determine how many **lines** are in /opt/bio-workshop/data/lamina.bed
    + use **wc** to determine how many **words** are in /opt/bio-workshop/data/lamina.bed
  

Less (is More)
--------------

To view a large file, use less:

.. code-block:: bash

    less /opt/bio-workshop/data/lamina.bed

You can forward-search in the file using "/"

You can backward-search in the file using "?"

You can see info about the file (including number of lines) using "ctrl+g"

You can exit **less** using "q"


Terminal History
----------------

Press the up arrow in the terminal.

Up and down arrows will allow you to scroll through your previous commands.

This is useful when running similar commands or when remembering what you have
done previously.

You can type the start of a command and then up-arrow and it will cycle
through commands that start with that prefix.


Tab-Completion
--------------

The shell (bash) when set up properly can give you a lot of help

Type the following where [TAB] means the Tab key on the keyboard:

.. code-block:: bash

    $ cd /opt/bio-w[TAB]

Then hit tab. And:

.. code-block:: bash

    $ ls /opt/bio-w[TAB]

This will work for any file path and for any programs.

.. code-block:: bash

    $ hea[TAB]

What happens if you do:

.. code-block:: bash

    $ he[TAB][TAB] 

or:

.. code-block:: bash


    $ heaaa[TAB][TAB] 


Directory Shortcuts
-------------------

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

Directory Shortcuts(2)
----------------------

We can go 2 directories up with:

.. code-block:: bash

    $ cd ../../

Here, we can remember that "." is the current directory and .. is one directory up.
What does this do:

.. code-block:: bash

    $ ls ./*

Directory Shortcuts(3)
----------------------

you can go to the last directory with:

.. code-block:: bash

    $ cd -

and switch back and forth by using that repeatedly.


make and remove directories
---------------------------

.. code-block:: bash

    mkdir ~/tmp # OK

    mkdir ~/tmp/asdf/asdf # ERROR

    mkdir -p ~/tmp/asdf/asdf # OK


What does -p do?

Remove directories:

.. code-block:: bash

   rm ~/tmp/asdf # ERROR

   rm -r ~/tmp/asdf/asdf # OK

What is -r ?

.. note ::

    be careful with `rm -r` and `rm -rf`

What does this do?

.. code-block:: bash

    rm -r /


moving/copying files
--------------------

mv [source] [dest]

.. code-block:: bash

    touch /tmp/asdf
    mv /tmp/asdf ~
    ls -lhtr ~/

moving/copying files(2)
-----------------------

in class excercise:


 1. make a directory `/tmp/moveable`
 2. move that directory to ~
 3. copy that directory to `/tmp/subdir/`


echo
----

echo is print:

.. code-block:: bash

    echo "hello world"

and you can use it to see **bash** variables:

.. code-block:: bash

    echo $HOME

    echo $HISTFILE

variables
---------

We will start covering programming in the next classes, but variables are a
key component of programming.

You can do:

.. code-block:: bash

    $important=/opt/bio-workshop/data/lamina.bed
    ls -lh $important


sudo
----

.. image:: http://imgs.xkcd.com/comics/sandwich.png

.. code-block:: bash

    apt-get install cowsay
    sudo apt-get install cowsay


other commands
--------------

excercise:

use `man` to determine the function of:

    + wget
    + uniq

How many records are present for each chromosome in
/opt/bio-workshop/data/lamina.bed (assume it is sorted by chromosome)?

gedit
-----

In order to edit files as you would using `notepad` or `word` in windows,
we will use the simple editor "gedit".

You can open gedit from the terminal using:

.. code-block:: bash

    $ gedit

This will open a new window with GUI controls. Use gedit to write/edit scripts for this class


Scripts
-------

A script is simply a series of commands that you save in a file. You will need to write
scripts to complete the homework.

Put this text:

.. code-block:: bash

    ls /opt/bio-workshop/

Into the file *`my-ls.sh`* by opening `gedit` pasting that text then `save as..` using the GUI controls

You can then run it as:

.. code-block:: bash

    bash my-ls.sh

And you should see the same output as if you ran `ls /opt/bio-workshop` directly.

Scripts
-------

Scripts will be more useful when you have a series of commands you want to run in series.

For example a pipeline where you:

 1. run quality control on some ChIP-seq reads 
 2. align reads to a reference genome
 3. find peaks (binding sites)
 4. annotate the binding sites.

In cases like that, a script will provide a record of what you have done.

Comments
--------

For the homework you will comment your scripts. 

Comments are not read by the shell, but they tell us (and you) what
you were trying to do. You can comment your code using the "#" symbol.

.. code-block:: bash
    
    # list all files in the /tmp/ directory ordered so that most recently
    # changed appear last
    $ ls -lhtr /tmp/

Pipes
-----

Since linux is made of small utilities, we often want to chain them
together. We will cover this in detail next class, but the idea
is that each program takes data, modifies it, and sends it to the next.

We can see lines 5-10 of a file with:

.. code-block:: bash

    head /opt/bio-workshop/data/lamina.bed | tail -n 5

Resources
---------

There is a nice summary of bash features here: http://digital-era.net/wp-content/uploads/2013/12/BASH-as-a-Modern-Programming-Language-Presentation-1.pdf

In Class Exercises
------------------


Place the answers to these in the bash script:


    1. write a bash script that you can run to list only the 2 most recently
       modified files in a given directory (using what you've learned in this class)
    2. make that script executable (use google to learn how to do this).

    3. With `head`, you can see the first line of a file with head -n1.
       How can you see all of a file *except* the first line.

    4. Without using your history, how few keystrokes can you use to run the following command (must work from any directory)?

        ls /opt/bio-workshop/data/lamina.bed

    5. How few keystrokes can you do #4 using your history?

    6. To learn about piping (|), use cowsay to:

       a. show your current working directory
       b. tell  you the number of lines in /opt/bio-workshop/data/lamina.bed
       c. tell you the most recently modified file (or directory) in $HOME
