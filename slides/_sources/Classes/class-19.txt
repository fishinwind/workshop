
.. include:: /_static/substitutions.txt

*************************
Class 19 : Applied Python 
*************************

:Class date: |c19-date|
:Last updated: |today|

Pybedtools
==========

The pybedtools documentation is very good: https://pythonhosted.org/pybedtools/

Citable from: http://www.ncbi.nlm.nih.gov/pubmed/21949271

.. code-block:: python

    for region in BedTool('my.bed'):
        if region.strand == "+":
            print region.chrom, region.start, region.end


.. raw:: pdf

    PageBreak
