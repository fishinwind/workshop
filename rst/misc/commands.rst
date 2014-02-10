***************************
Reference: list of commands
***************************

cd
==
change directories:

.. code-block:: bash

    $ cd /tmp/
    $ cd ~ # chage to home directory
    $ cd /opt/bio-workshop/
    $ cd - # change to previous directory

cp
==
copy files and directories:

.. code-block:: bash

    $ touch /tmp/asdf
    $ cp /tmp/asdf ~ # copy to home

must use -r for directories:

.. code-block:: bash

   $ mkdir /tmp/adir
   $ cp -r /tmp/adir ~/

ctrl+c
======
interrupt a running process:

.. code-block:: bash

    $ head
    <ctrl+c>

cut
===
extract columns from a file:

.. code-block:: bash

    $ cut -f 1-3 /opt/bio-workshop/data/lamina.bed
    $ cut -f 1,3 /opt/bio-workshop/data/lamina.bed
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed

    # use comma as delimiter instead of default tab
    $ cut -f 1-3 -d, /path/to/some.csv

    # keep all columns after the 1st:
    $ cut -f 2- /opt/bio-workshop/data/lamina.bed

echo
====
print some text:

.. code-block:: bash

    $ echo "hello world" | cowsay

head
====
show the start of a file:

.. code-block:: bash

    $ head /opt/bio-workshop/data/lamina.bed

    # show the first 4 lines
    $ head -n 4 /opt/bio-workshop/data/lamina.bed

grep
====
To find any instance of *chr5* in the lamina.bed file

.. code-block:: bash

    # grep [pattern] [filename]
    $ grep chr5 /opt/bio-workshop/data/lamina.bed | head

To find all lines that start with a number sign:

.. code-block:: bash

    # The caret (^) matches the beginning of the line
    # FYI dollar sign ($) matches the end
    $ grep '^#' /opt/bio-workshop/data/lamina.bed

To find any line that *does not* start with "chr":

.. code-block:: bash

    # the -v flag inverts the match (grep "not" [pattern])
    $ grep -v '^chr' /opt/bio-workshop/data/lamina.bed

Find exact matches that are split on words with the ``-w`` flag:

.. code-block:: bash

    $ grep -w chr1 /opt/bio-workshop/data/lamina.bed | cut -f1 | uniq

less
====
page through a file:

.. code-block:: bash

    $ less /opt/bio-workshop/data/lamina.bed

use "/", "?" to search forward, backard. 'q' to exit.

use '[space]' to go page by page.

ls
==
list files and directories:

.. code-block:: bash

    $ ls /tmp/

    # show current directory
    $ ls

    # show current directory (2)
    $ ls . 

    # list files with most recently modified last
    $ ls -lhtr

    # list files in temp ordered by modification date
    $ ls -lhtr /tmp/ 

man
===
show the manual entry for a command:

.. code-block:: bash

    $ man head

mkdir
=====
make a directory:

.. code-block:: bash

    $ mkdir /tmp/aa

make sub-directories, too:

.. code-block:: bash

    $ mkdir -p /tmp/aaa/bbb/

mv
==
move a file or directory:

.. code-block:: bash

    $ touch /tmp/aa
    $ mv /tmp/aa /tmp/bb

rm
==
remove a file or directory:

.. code-block:: bash

    $ touch /tmp/asdf
    $ rm /tmp/asdf

    # use -r to remove directory
    $ mkdir /tmp/asdf
    $ rm -r /tmp/asdf

sort
====
sort a file by selected columns:
    
.. code-block:: bash

    $ sort -k1n /opt/bio-workshop/data/lamina.bed

sort a BED file by chromosome (1st column) as character and then by start
(2nd column) as number:

.. code-block:: bash

    $ sort -k1,1 -k2,2n /opt/bio-workshop/data/lamina.bed

sort by 4th column as a general number, including scientific notation
showing largest numbers first:

.. code-block:: bash

    $ sort -k4,4rg /opt/bio-workshop/data/lamina.bed | head

use literal tab ('\\t') as the delimiter (default is whitespace):

.. code-block:: bash

    $ sort -t$'\t' -k4,4rg /opt/bio-workshop/data/lamina.bed | head

Sometimes we want to get uniq entries with sort -u:

.. code-block:: bash

    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | sort -u

will print out the uniq chromsomes represent in the BED file.

tail
====
show the end of a file:

.. code-block:: bash

    $ tail /opt/bio-workshop/data/lamina.bed
    # show the last 4 lines
    $ tail -n 4 /opt/bio-workshop/data/lamina.bed

tar
===
create or untar a .tar.gz file:

.. code-block:: bash

    # -c create -z compress (.gz) -v verbose -f the name 
    $ tar -czvf some.tar.gz /tmp/*

    # -x untar 
    $ tar -zxvf some.tar.gz

uniq
====
show or count unique or non-unique entries in a file:

.. code-block:: bash

    # count number of times each chromosome appears.
    $ cut -f 1 /opt/bio-workshop/data/lamina.bed | uniq -c

    # get non unique entries
    $ cut -f 2 /opt/bio-workshop/data/lamina.bed | uniq -d

.. important::
    
    ``uniq`` assumes that the file is ``sort``-ed first! To test this, run
    ``uniq`` on an unsorted file. What happens?

wget
====
get a file from the web:

.. code-block:: bash

    $ wget http://ucd-bioworkshop.github.io/_downloads/states.tab

zless
=====
like less, but for compressed files:

.. code-block:: bash

    $ zless /opt/bio-workshop/data/t_R1.fastq.gz

Redirection (>> and >)
======================
send output to a file:

.. code-block:: bash

    # start a new file
    $ echo "hello" > file.txt 

    # overwite that file
    $ echo "hello!" > file.txt 

    # append to the file
    $ echo "world" >> file.txt 

.. raw:: pdf

    PageBreak
