***************
Software for VM
***************

VM Installation
===============
- The current VM is ~4.5 Gb in size, which means it can not be copied to a
  FAT-formatted USB stick. Reformat the stick with exFAT and all is well
- VirtualBox gives funny errors on some personal laptops. Have folks check
  whether Virtualization is an option in the BIOS, and if so, enable it.

Software for the VM
===================

General
-------
- gedit
- Rstudio

Bioinformatic tools
-------------------
- bwa: https://github.com/lh3/bwa
- samtools: https://github.com/lh3/samtools
- bioawk: https://github.com/lh3/bioawk
- tophat2:  http://tophat.cbcb.umd.edu/downloads/tophat-2.0.10.tar.gz
- bedtools: https://github.com/arq5x/bedtools
- fastqc: (zip and .jar) http://www.bioinformatics.babraham.ac.uk/projects/download.html#fastqc

R packages
----------
 - ggplot2
 - limma
 - biocondcutor

Python 2.7+ packages
--------------------
 - numpy (apt)
 - scipy (apt)
 - matplotlib (apt)
 - pandas (pip)
 - toolshed (pip)

VM installs via apt-get
-----------------------
 - build-essential
 - zlib1g-dev
 - r-base-core, r-base
 - bash-completion
 - git
 - svn
 - mercurial 

