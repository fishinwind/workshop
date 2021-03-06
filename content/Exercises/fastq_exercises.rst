
.. include:: /_static/substitutions.txt 

.. _exercises:

*************
  FASTQ/A Exercises
*************

:Last updated: |today|

Exercises using unix tools to manipulate FASTQ/A files 

.. _fastq-exercises:

FASTQ exercises
==============

Download a fastq example file :download:`[SP1.fq] <../Miscellaneous/data/SP1.fq>`

The FASTQ file `format <https://en.wikipedia.org/wiki/FASTQ_format>`_ contains 4 lines per read sequence. 

.. code-block:: bash

    $ head -n 4 SP1.fq
    @cluster_2:UMI_ATTCCG
    TTTCCGGGGCACATAATCTTCAGCCGGGCGC
    +
    9C;=;=<9@4868>9:67AA<9>65<=>591


#. Line 1 contains a sequence identifier that begins with an `@`
#. Line 2 contains the read sequence (A,T,C,G,N) 
#. Line 3 is a comment line, often unused and only contains a `+` 
#. Line 4 is the Phred quality score for each sequence encoded in ASCII format 


First use ``wc`` and ``awk`` to determine the number of *sequences* in the fastq 

.. code-block:: bash
    
    $ wc -l SP1.fq  # total number of lines
    $ wc -l SP1.fq | awk '{print $1 / 4}' # total number of fastq records


A common mistake is to use ``grep`` to pattern match the ``@`` in the
sequence identifier. Why doesn't this work?

.. code-block:: bash
    
    $ wc -l SP1.fq | awk '{print $1 / 4}'
    250
    $ grep -c "@" SP1.fq
    459


Next, extract out the read *sequences* from the fastq. This is a bit
complicated as we need to only pull out the second line from each record. 

One approach to this problem is to use the ``%`` `modulo
operator <https://en.wikipedia.org/wiki/Modulo_operation>`_, which returns
the remainder after division of two integers. For example using ``awk``:   

.. code-block:: bash

    $ awk 'BEGIN { {print 4 % 2}}'
    $ awk 'BEGIN { {print 4 % 3}}'
    $ awk 'BEGIN { {print 5 % 3}}'
    $ awk 'BEGIN { {print 1 % 4}}'

In ``awk`` there is a special variable ``NR`` which is equal to the
current line number.


.. code-block:: bash

    $ awk '{print NR}' SP1.fq


We can use the modulo operator with ``NR`` to only capture specific
records from the fastq

.. code-block:: bash

    $ awk '{print NR % 4, $0 }' SP1.fq | head # note output in first column
    1 @cluster_2:UMI_ATTCCG
    2 TTTCCGGGGCACATAATCTTCAGCCGGGCGC
    3 +
    0 9C;=;=<9@4868>9:67AA<9>65<=>591
    1 @cluster_8:UMI_CTTTGA
    2 TATCCTTGCAATACTCTCCGAACGGGAGAGC
    3 +
    0 1/04.72,(003,-2-22+00-12./.-.4-
    1 @cluster_12:UMI_GGTCAA
    2 GCAGTTTAAGATCATTTTATTGAAGAGCAAG


    $ awk 'NR % 4 == 1' SP1.fq | head  # get header line
    $ awk 'NR % 4 == 2' SP1.fq | head  # get sequence line
    $ awk 'NR % 4 == 3' SP1.fq | head  # get comment line
    $ awk 'NR % 4 == 0' SP1.fq | head  # get quality line


Now that we can isolate only the sequences let's use `sort` and `uniq` to
find common sequences. 

.. code-block:: bash

    awk 'NR % 4 == 2' SP1.fq \  # get sequences
      | sort \  # sort the sequences 
      | uniq -c \ # report number of occurances of each unique sequence
      | head 

    awk 'NR % 4 == 2' SP1.fq \  #show most frequent sequences
      | sort \
      | uniq -c \
      | sort -k1,1nr \ # reverse sort to rank by most common sequence
      | head 

Next, using ``cut``, ``sort`` and ``uniq`` lets look at the most common
sequences found in the first 10 base pairs of the read.

.. code-block:: bash

    awk 'NR % 4 == 2' SP1.fq \  # get sequences
    | cut -c 1-10 \ # cut by character position 1-10
    | sort \
    | uniq - c \
    | sort -k1,1nr \
    | head 

``awk`` can also be used directly to extract out sub-strings using the ``substr()``
function

.. code-block:: bash

    # substr($field_index, start, end)
    $ awk 'NR % 4 == 1 {print substr($0, 1, 10)}' SP1.fq \
      | head -5 

Let's extract out the Unique Molecular Identifyier (UMI)
tag from the read identifier using ``cut`` 


.. code-block:: bash

    $ awk 'NR % 4 == 1' SP1.fq | head -5 # get ids
    @cluster_2:UMI_ATTCCG
    @cluster_8:UMI_CTTTGA
    @cluster_12:UMI_GGTCAA
    @cluster_21:UMI_AGAACA
    @cluster_29:UMI_GCAGGA

    $ awk 'NR % 4 == 1' SP1.fq \
      | cut -d ":" -f 2 \
      | head -5 


Lastly, let's write a short awk program that uses the ``substr()``
function to trim the first 10 bases off of the read and return the
sequences in FASTQ. Note that the length of the quality line needs to
match the length of the sequence line. 

.. code-block:: bash

    $ awk 'NR % 4 == 1 {print $0}; # print id line
      NR % 4 == 2 {print substr($0, 10, length($0))}; #extract from 10 to the end 
      NR % 4 == 3 {print $0} ; # print comment line
      NR % 4 == 0 {print substr($0, 10, length($0))}' \ # extract from pos 10 to the end 
      SP1.fq \
      | head 

    @cluster_2:UMI_ATTCCG
    CACATAATCTTCAGCCGGGCGC
    +
    4868>9:67AA<9>65<=>591
    @cluster_8:UMI_CTTTGA
    AATACTCTCCGAACGGGAGAGC
    +
    003,-2-22+00-12./.-.4-
    @cluster_12:UMI_GGTCAA
    GATCATTTTATTGAAGAGCAAG


.. _fasta-exercises:

FASTA exercises
==============

FASTA format is a more compact sequence format than FASTQ

.. code-block:: bash

    >sequence_identifier_1
    ATCGTCGATGCTAGTCGA
    >sequence_identifier_2
    AGCTAGCTAGCTAGC

Using just ``awk`` we can easily convert from FASTQ to FASTA

.. code-block:: bash

    $ awk 'NR % 4 == 1 {print ">"$1}; 
      NR % 4 == 2 {print}' SP1.fq \
      | head 

    $ awk 'NR % 4 == 1 {print ">"substr($0, 2, length($0))}; 
      NR % 4 == 2 {print}' SP1.fq \
      | head  # same as above but remove ugly @ character

    $ awk 'NR % 4 == 1 {print ">"substr($0, 2, length($0))}; 
      NR % 4 == 2 {print}' SP1.fq \
      > example.fa 

To count the number of sequences in the fasta file you can use either
``wc`` with ``awk`` as shown above or ``grep`` as follows:

.. code-block:: bash
 
    grep -c ">" example.fa 

On your own
==========

#. Use ``awk`` to extract out sequences and find out the lengths of each
   read. Are they all the same length? 
#. Use ``grep`` and/or ``awk`` to find sequences with more than 4 A's in a row. How many
   sequences did you find? Can you also extract out their names (hint see
   man grep)? Do any UMI's have 4 A's in a row?
#. What is the most common hexamer found at the end of the reads? 
