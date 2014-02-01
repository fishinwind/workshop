Week 3 / Class 1 : Working in a cluster environment
===================================================

Cluster access
--------------
We have set up accounts for the class on our departmental cluster. We will
set up your accounts at the end of class and reset your passwords.

.. code-block:: bash
    
    $ ssh -X username@amc-tesla.ucdenver.pvt

Cluster etiquette
-----------------
There are some specific rules you need to know when you're operating in a
cluster environment.

- The cluster is run by a central computer called the *head node*. This is
  the computer that you log into (amc-tesla). **DO NOT** run jobs on the
  head node. The head node is essentially the brains of the cluster, and
  it can easily be overextended.

The queueing system
-------------------
The cluster uses a queueing system that will run jobs that you submit to
it. You can write a small test script to see how the system works. First,
write this into a run.sh file:

.. code-block:: bash

    # /usr/bin/env bash

    #BSUB -J sleeper
    #BSUB -e %J.err
    #BUSB -o %J.out

    sleep(20)

The `#BSUB` lines are comments, but are read by the submission program to
identify features associdated with your jab. The `-J` flag sets the jobs
name. The `-e` and `-o` options set the filenames for the output from the
job printed to stdout and stderr

Now you can submit the script to the queuing system. As soon as you submit
it, you can check on its progress.

.. code-block:: bash

    $ bsub < run.sh
    $ bjobs

Killing jobs
------------
Sometimes you need to kill your jobs. You can kill specific jobs using
their job ID numbers, obtained from checking `bjobs`:

.. code-block:: bash

    $ bkill <jobid>

You can also kill **all** of your jobs at once. **Use this sparingly**.
If you have long running jobs in addition to jobs you just submitted,
you'll wipe out everything:

.. code-block:: bash

    $ bkill 0


