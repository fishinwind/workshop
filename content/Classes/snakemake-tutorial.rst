
.. include:: /_static/substitutions.txt

=============================
Class 18 : Reproducible Research using Snakemake 
=============================

:Last updated: |today|

Goals
-----
#. Introduce reproducible research methods using pipelines and 
   `Snakemake <https://snakemake.readthedocs.io/en/latest/index.html>`_
#. Generate Snakemake pipeline for Chip-Seq analysis

Overview
--------

#. Genomics analysis projects often generate hundreds or thousands of files.
   Keeping track of all of the scripts used to generate a set of finalized
   reports is an error-prone and difficult task.

#. The naive approach to organizing a project is to make a ``bash`` script
   that executes each step of the analysis from mapping to plotting. 

#. Organizing a project in this manner is a good first step to producing
   reproducible analyses, however there are few common problems with this
   approach. 

For example:

#. After plotting all of your results you realize that you need to rerun
   the peak calling step, but in order to do that you need to rerun the
   entire pipeline.

#. The script might grow to be 1,000s of lines of code making debugging
   very tedious as mistakes will be hard to spot. 

Modularizing code
-----------------
The next approach might be to split each important set of steps into
individual scripts. i.e:

#. ``run_mapping.sh``

#. ``call_peaks.sh``

#. ``plot_data.R``

By modularizing the code, the workflow is easier to modify and it
becomes easier to reuse the code. However, with this setup a few
additional problems can occur. 

Potential Problems
------------------

#. You rerun the ``run_mapping.sh`` script, but forget to rerun
   ``call_peaks.sh`` and ``plot_data.R``. Two weeks later you discover your
   mistake...  

#. You edit and rerun ``call_peaks.sh`` but an error occurs, and now you are unsure
   which of the hundreds of files you generated are new or old.  

For these reasons, it is a very good practice to use a workflow management
system to help with managing computational projects.

Installation
===========

.. code-block:: bash
    
    pip3 install snakemake

Snakemake
=========

`Snakemake <https://snakemake.readthedocs.io/en/latest/index.html>`_ is a
worflow manangement system that is based upon GNU Make and written in python.

The central idea of a snakemake workflow is that we define a set of rules
that specify how to generate output files from input files. Snakemake uses
these rules to track which files need to be generated, and if the files are missing
will automagically regenerate files as necessary.

Hello World
-----------
Let's start with a simple example:

.. code-block:: bash

    # generate an example file
    $ echo "world" > world.txt


Copy the following code into a file named ``Snakefile`` (no extension).

.. code-block:: basemake

    rule hello_world:
      input: "world.txt"
      output: "hello_world.txt"
      shell:
        " echo 'Hello' |  cat - {input} > {output} "

Run Snakemake
=============
Now let's use Snakemake to execute our rule. 

.. code-block:: bash
    
    snakemake -npr hello_world.txt 
    # -npr tells snakemake to do a dry run and print the expected commands

Rerun Snakemake
===============
As we can see snakemake has filled in the ``{input}`` and ``{output}`` variables
for our rule. Next go ahead and exectute this simple rule (remove
``-npr``). 

Now here's where snakemake get's useful. Go ahead and change the contents
of the ``world.txt`` file. Now when we reexecute Snakemake, it knows that
``world.txt`` has been edited, and that ``hello_world.txt``, which depends
on ``world.txt`` also needs to be regenerated


.. code-block:: bash
    
    snakemake -npr hello_world.txt 

Snakemake will automatically reexecute a rule if the input files have been
modified. If we have many rules, then any intermediate files will also be
regenerated. For example, add another rule that now edits
``hello_world.txt``:

Hello Universe
==============
.. code-block:: basemake

    rule hello_universe:
      input: "hello_world.txt"
      output: "hello_universe.txt"
      shell:
        "sed 's/world/universe/' {input} > {output} "

.. code-block:: bash
    
    snakemake -npr hello_universe.txt 

Rule All
========

#. Instead of defining our target files on the command line we can specify
   them directly in the ``Snakefile``. 
#. Rule ``all`` is a psuedo-rule that simply tells snakemake what final files to generate. 
#. By default snakemake determines which rule to execute first in a top-down manner. 
#. Therefore snakemake will first search for ``hello_universe.txt``, then will find the rule
   that generates it (``rule hello_universe``), then execute any additional rule
   dependencies, until all target files have been generated. 

Rule All Example
================

.. code-block:: basemake

    rule all:
      input: "hello_universe.txt"
    
    rule hello_universe:
      input: "hello_world.txt"
      output: "hello_universe.txt"
      shell:
        "sed 's/world/universe/' {input} > {output}
    
    rule hello_world:
      input: "world.txt"
      output: "hello_world.txt"
      shell:
        "echo 'Hello' | cat - {input} > {output}"

Generalizing rules with wildcards
=================================

Snakemake supports using wildcards to define the ``input`` and ``output``
files. 

.. code-block:: basemake

    PLACES = ["hello_universe.txt",
              "hello_colorado.txt",
              "hello_123.txt"]
    rule all:
        input:
          PLACES #list of files is substituted
    rule hello_universe:
        input: "hello_world.txt"
        output: "hello_{place}.txt" 
        shell:
          "sed 's/world/{wildcards.place}/' {input} > {output}"
          
    rule hello_world:
      input: "world.txt"
      output: "hello_world.txt"
      shell:
        "echo 'Hello' | cat - {input} > {output} "

Expand Function
===============
Snakemake has a helpful python function called ``expand`` that simplifies
the previous code to:

.. code-block:: basemake

    PLACES = ["universe",
              "colorado",
              "123"]

    rule all:
        input:
          expand("hello_{place}.txt", place=PLACES) #list of files is substituted 

Debugging Targets
=================

To check that your variable substitutions are correct you can use the
python ``print()`` function, which will print the substitued variables to
standard out, or alternatively look at the snakemake output (``-npr``)

.. code-block:: basemake

    PLACES = ["universe",
              "colorado",
              "123"]

    print(expand("hello_{place}.txt", place=PLACES)) 


Snakemake pipeline for Chip-Seq analysis
----------------------------------------

So far our current pipeline isn't very useful. Next we'll write a short
pipeline to conduct basic Chip-Seq analysis. In this case we will use
H3k4me3 Chip-Seq data from hESC cells generated by the encode project,
with the goal of examining the distribution of H3k4me3 around human
Transcriptional Start Sites (TSS). 

Important steps for the pipeline: 

#. Map the data to the genome using ``bowtie2`` (FASTA -> BAM)

#. Sort and Index `BAM` (BAM -> Sorted BAM)

#. Get coverage using ``bedtools genomecov`` (BAM -> bedGraph)

#. Plot coverage over TSS using bedtools/R` (bedGraph -> pdf)

Some of the steps in the pipeline were covered in the Bedtools vignette class.

Get necessary software
======================

First lets make sure that we have all of the necessary software installed. 

.. code-block:: bash

    brew install bowtie2 samtools bedtools

Make a new file named ``Snakefile`` in the ``plot_coverage`` directory. 

Define global variables
======================
Let's start our pipeline by defining a few important variables
for our analysis. Place these at the top of your Snakefile. 

.. code-block:: basemake

    FASTA = "dbases/chr22.fa.gz" # our genome fasta file (hg19)
    CHROMS = "dbases/chr22_length.txt" # our genome file necessary for bedtools
    GENES = "dbases/knownGene_chr22.bed" # gene annotations for hg19


It is convienent to organize the pipeline in a bottom-to-top fashion,
whereby the first steps executed are at the bottom. 

For this analysis we will first align the Chip-Seq data using ``bowtie2``.
First we'll generate the index necessary for alignment, then perform the
alignment. 

Indexing 
========

.. code-block:: basemake

   rule bowtie_index:
     input: FASTA 
     output: "dbases/bowtie_idx/chr22.1.bt2"
     params:
       output_name = "dbases/bowtie_idx/chr22"
     shell:
       """
       gunzip -c {input} > {input}.tmp
       bowtie2-build {input}.tmp {params.output_name} 
       """

Alignment
=========
.. code-block:: basemake

   rule bowtie_mapping:
     input: 
       fq = "raw_data/{chip}.fastq.gz", 
       idx = "dbases/bowtie_idx/chr22.1.bt2"
     output: "bowtie/{chip}.bam"
     params: idx = "dbases/bowtie_idx/chr22", 
     shell:
       """
       bowtie2 \
         -x {params.idx} \
         -U {input.fq} \
         -S {output}.tmp 

       samtools sort {output}.tmp > {output}
       samtools index {output}
       """

Key Concepts
============

#. multiple ``inputs`` or ``outputs`` can be defined and called by name
   i.e. (``{input.fq}``). Each input must be separated by a comma.

#. Arbitrary additional arguments can be specified using the ``params``
   option. This is useful for customizing the command line arguments
   with variables, but not requiring input or output dependencies. 

#. We have a wildcard (``chip``) that takes the place of the fastq name.
   This wildcard will allow us to generalize our pipeline. 

#. Global variables, such as ``FASTA`` can be used as an input. Global
   variables can also be used in the ``shell:`` command but need to be
   bracketed (``{FASTA}``). 

#. Snakemake will automatically generate directories that are listed in
   the output files, if they do not exist. 

Check Alignment Step
====================

Check that our snakemake pipeline is working:

.. code-block:: bash
    
    snakemake -npr bowtie/H1_h3k4me3_chr22.bam

Snakemake will pattern match our requested file
(``bowtie/H1_h3k4me3_chr22.bam``) and determine that it can make this file
by substituting ``H1_h3k4me_chr22`` for the ``chip`` variable. 

Coverage
========

Next calculate alignment coverage across the genome with bedtools.

.. code-block:: basemake

    rule make_bedgraphs:
      input:
        "bowtie/{chip}.bam"
      output:
        "bowtie/{chip}.bedgraph"
      shell:
        """
        bedtools genomecov \
          -ibam {input} \
          -bg \
          -g {CHROMS} > {output}
        """


Also note that we used the global variable ``CHROMS`` directly in the
shell command. The braces are needed in the shell command so that
snakemake can recognize it, but are not needed when calling the variable
outside of the shell commands. 

Check Coverage
==============

Let's check that our pipeline is working by trying to generate the
bedgraph file. 

.. code-block:: bash
    
        snakemake -npr bowtie/H1_h3k4me3_chr22.bedgraph

Make TSS windows
================
Now that we have our bedgraph we next need to define a set of intervals to
surround the TSS to calculate H3k4me3 coverage. To do this we will call an
external bash script that uses awk and bedtools to generate a bed file of
windows around each TSS. 

.. code-block:: basemake

    rule prepare_windows:
      input:
        genes = GENES,
        chroms = CHROMS,
      output:
        "bedfiles/tss_windows.bed"
      shell:
        "src/get_tss_windows.sh {input.genes} {input.chroms} {output}" 
    
These rules take the reference gene annotations and extracts out the TSS,
then make a set of windows surrounding the TSS.

Other Concepts
==============
#. Instead of writing an external script, you can also include each step
   individually. However, these steps are very verbose, so they have been
   moved to a seperate script.  

#. Snakemake doesn't care what commands we run, only that the input files are
   present before running and that the output files are generated after
   executing the commands. 

Let's check that our pipeline is correct by trying to generate the
output from the ``prepare_windows`` rule. 

.. code-block:: bash
    
        snakemake -npr "bedfiles/tss_windows.bed"


Calculate coverage over TSSs
============================

Next let's use bedtools map to calculate the mean H3K4me3 signal over
each TSS region

.. code-block:: basemake

    rule get_coverage:
      input: 
        windows = "bedfiles/tss_windows.bed",
        bedgraph = "bedgraphs/{chip}.bedgraph"
      output: 
        coverage = "bedfiles/{chip}/coverage.bed"
      shell:
        """
        bedtools map -a {input.windows} \
          -b {input.bedgraph} \
          -c 4 -o mean -null 0 > {output.coverage}
        """

Check bedtools map step
=======================
Again check the commands:

.. code-block:: bash
    
    snakemake -npr "bedfiles/H1_h3k4me3/coverage.bed"

Plot the Data
=============

Lastly, let's group and plot the data using the `bin/plot_data.R` script. 

.. code-block:: basemake

    rule plot_coverage:
      input:
        coverage = "bedfiles/{chip}/coverage.bed"
      output:
        "plots/{chip}.pdf"
      shell:
        " bin/plot_data.R {input} {output} "

Again check the commands:

.. code-block:: bash
    
    snakemake -npr "plots/H1_h3k4me3_chr22.pdf"


Using Rule All
==============
To clean up the pipeline let's define our target files in the Snakefile
itself. 

.. code-block:: basemake
    
    CHIP = ["H1_h3k4me3_chr22"]
    
    rule all:
      input: expand("plots/{chip}.pdf", chip = CHIP)

.. code-block:: bash
    
    snakemake -npr 

Putting it all together
=======================

.. code-block:: basemake

    FASTA = "dbases/chr22.fa" # our genome fasta file (hg19)
    CHROMS = "dbases/chr22_length.txt" # our genome file for bedtools
    GENES = "dbases/knownGene_chr22.bed" # gene annotations for hg19
    
    CHIP = ['H1_h3k4me3_chr22'] 
    
    rule all:
      input: expand("bedfiles/{chip}/coverage.bed", chip=CHIP)
    
    rule plot_coverage:
      input:
        coverage = "bedfiles/{chip}/coverage.bed"
      output:
        "plots/{chip}.pdf"
      shell:
        " src/plot_data.R {input} {output}"
    
    rule get_coverage:
      input: 
        windows = "bedfiles/tss_windows.bed",
        bedgraph = "bedgraphs/{chip}.bedgraph"
      output: 
        coverage = "bedfiles/{chip}/coverage.bed"
      shell:
        """
        bedtools map -a {input.windows} \
          -b {input.bedgraph} \
          -c 4 -o mean -null 0 > {output.coverage}
        """
    
    rule prepare_windows:
      input:
        genes = GENES,
        chroms = CHROMS,
      output:
        "bedfiles/tss_windows.bed"
      shell:
        "src/get_tss_windows.sh {input.genes} {input.chroms} {output}" 
        
    rule make_bedgraphs:
      input:
        "bowtie/{chip}.bam"
      output:
        "bedgraphs/{chip}.bedgraph"
      shell:
        """
        bedtools genomecov \
          -ibam {input} \
          -bg \
          -g {CHROMS} > {output}
        """
    
    rule bowtie_mapping:
      input: 
        fq = "raw_data/{chip}.fastq.gz", 
        idx = "dbases/bowtie_idx/chr22.1.bt2"
      output:
        "bowtie/{chip}.bam"
      params:
        idx = "dbases/bowtie_idx/chr22", 
      shell:
        """
        bowtie2 \
          -x {params.idx} \
          -U {input.fq} \
          -S {output}.tmp
    
        samtools sort {output}.tmp > {output}
        samtools index {output}
        """
    
    rule bowtie_index:
      input: 
        FASTA 
      output:
        "dbases/bowtie_idx/chr22.1.bt2"
      params:
        output_name = "dbases/bowtie_idx/chr22"
      shell:
        """
        bowtie2-build {input} {params.output_name} 
        """

Generalizing the pipeline
=========================

The real power of a pipeline is the ability to apply it to many datasets.
By adding additional elements to our ``CHIP`` variable we can now run our
pipeline on multiple Chip-Seq samples

.. code-block:: basemake
    
    CHIP = ["H1_h3k4me3_chr22", "H1_RNApolII_chr22"]

Now our pipeline will run for both Chip-Seq datasets with only 1 filename
added!

Glob wildcards
==============

#. What if you have a directory with 100s of fastq files from Chip-Seq
   experiments? You could type out a long list of filenames, or ...

#. Use the ``glob_wildcards`` function in Snakemake. This function will build a
   list of sample names based on the files in the listed directory.  

.. code-block:: basemake

    CHIP, = glob_wildcards("raw_data/{chip}.fastq.gz")

Run in Parallel
==============

Lastly, snakemake can run jobs in parallel. If your computer has more than
1 core you can execute more than 1 job at the same time. Snakemake handles
the scheduling, all you need to do is tell it the number of cores you want
to use.

.. code-block:: bash
    
    snakemake --cores 2 

Troubleshooting
===============
https://snakemake.readthedocs.io/en/latest/project_info/faq.html

http://stackoverflow.com/questions/tagged/snakemake

https://groups.google.com/forum/#!forum/snakemake

Helpful guidelines
==================

#. Indentation is important, use two or four spaces for each indentation.
#. Define your target (final output) files in rule all
#. Use unique extensions or directories for each rule to avoid wildcard
   collisions

Alternative tutorials
=====================
https://snakemake.readthedocs.io/en/latest/index.html

https://github.com/leipzig/SandwichesWithSnakemake

Exercises
=========

#. Make a three rule snakemake pipeline using basic unix tools. Include at
   least one wildcard variable. 

#. Download new Chip-Seq fastqs into the following directory from `Encode <https://www.encodeproject.org/matrix/?type=Experiment>`_. 
   Modify the snakemake pipeline to run the analysis on all the new fastqs. What do
   you notice about the TSS coverage for each Chip-Seq expt? 

#. Add additional rules to our pipeline to plot Chip-Seq coverage at the
   3' transcriptional stop site.

