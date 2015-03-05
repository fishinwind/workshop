Python Debugging
================
Many people are struggling with simple debugging issues. Stratgies:

Use :py:func:`print` copiously. For example, Print out a nested
data structure as it is being created

.. code-block:: python
    :emphasize-lines: 10

    from collections import defaultdict

    nested = defaultdict(list)

    for chrom, start, end in parse_bed(bedfile):

        coords = (start, end)
        nested[chrom].append(coords)

        print nested

.. nextslide::
    :increment:

Use the Python Debugger: :py:mod:`pdb`.

.. code-block:: python
    :emphasize-lines: 1,5

    import pdb

    for idx in range(100):
        ...
        pdb.set_trace()

    # keys:
    # c(ontinue) : move through the current context (i.e. loop)
    # <ctrl>-c   : exit the debugger
    (Pdb) idx
    0
    (Pdb) c
    (Pdb) idx
    1


