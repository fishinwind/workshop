***********
Class Notes
***********

Notes for Spring 2015
=====================

- Merge awk class into Linux command line classes, cover grep much earlier,
  and trim this to squeeze in more Python

- Students struggled with problem set 3 (the first python problem set).
  we failed to go into detail on nested data structures (e.g. dict of
  lists of tuples). we also probably should have gone over multi-line
  data formats a bit more prior to the problem set.

- Next year, should focus less on trivial examples (e.g. fruit examples) and
  spend more time going over biologic data.

#. create nested data types during loops
#. logic for determining state during looping
#. functions to act on data structures once created:

- min(), max()
- Counter.most_common()

- In class 10: Dont included the extra lines in `sample-seq-info.csv` it
  was too confusing to use the skip_while arg to toolshed.reader to pass a
  function

Name expansion needs to be explained a bit better. one slide should have:

.. code-block:: python

    # these are all the same ...

    # explicit but cluttered
    for tuple in zip(list1, list2):
        item1 = tuple[0]
        item2 = tuple[1]
        ...

    # better, parens not needed
    for (item1, item2) in zip(list1, list2):
        ...

    # best
    for item1, item2 in zip(list1, list2):
        ...

    # another tuple expansion 
    for idx, item in enumerate(list):
        ...

- should include ``BEDTools`` examples from
  https://github.com/arq5x/bedtools-protocols

