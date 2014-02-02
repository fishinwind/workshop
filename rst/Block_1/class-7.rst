******************************************
Class 7 : Working in a cluster environment
******************************************

Goals
=====
#. Login to cluster
#. Learn about cluster-specific commands
#. Queueing system basics

Cluster access
==============
We have set up accounts for the class on our departmental cluster. We will
set up your accounts at the end of class and reset your passwords::

    # the -X flag starts an X11 connection 
    $ ssh -X username@amc-tesla.ucdenver.pvt

Cluster etiquette
=================
There are some specific rules you need to know when you're operating in a
cluster environment.

- The cluster is run by a central computer called the *head node*. This is
  the computer that you log into (amc-tesla). **DO NOT** run jobs on the
  head node. The head node is the brains of the cluster and
  it can easily be overextended.

Example commands on the cluster
===============================
Find the size of the file system::

    $ df -h

Find how much space you have allocated::

    $ quota -h

The queueing system
===================
First you will grab a single CPU from the queueing system so that you can play
around without affecting the head node. We use ``qlogin`` for this::

    jhessel@amc-tesla ~
    $ qlogin 

    Job <492536> is submitted to queue <interactive>.
    <<ssh X11 forwarding job>>
    <<Waiting for dispatch ...>>
    <<Starting on compute00>>

    jhessel@compute00 ~
    $ 

.. note:: 

    The host in the prompt changed from ``amc-tesla`` to ``compute00``.
    
You can now execute long-running processes without worry of affecting the
cluster. Type ``exit`` to return back to your head node login.

The queueing system (2)
-----------------------
The cluster uses a queueing system that will run jobs that you submit to
it. You can write a small test script to see how the system works. First,
write this into a run.sh file::

    # /usr/bin/env bash

    #BSUB -J sleeper
    #BSUB -e %J.err
    #BUSB -o %J.out

    sleep(20)

The `#BSUB` lines are comments, but are read by the ``bsub`` program to
identify features associated with your job. 

    - ``-J`` sets the job's name
    - ``%J`` is a unique job ID that is set when you run the job.
    - ``-e`` and ``-o`` set the filenames for stderr and stdout from the job

The queueing system (3)
-----------------------
Now you can submit the script to the queuing system. As soon as you submit
it, you can check on its progress::

    $ bsub < run.sh
    $ bjobs

After the job finishes, you should see two new files that end
`.out` and `.err`; these stdout and stderr from the running job.
Look at the contents of those files so you know what is in
each one.

Killing jobs
============
Sometimes you need to kill your jobs. You can kill specific jobs using
their job ID numbers, obtained from checking `bjobs`::

    $ bkill <jobid> 

You can also kill **all** of your jobs at once::

    $ bkill 0 

.. warning::

    ``bkill 0`` is dangerous – it will wipe out all of your jobs. If
    you have long-running jobs that you forgot about, you will kill them
    too if you are not careful!

Other cluster-specific commands
===============================
.. code-block:: bash

    $ bhosts     # hosts in the cluster
    $ man bhosts # bsub man page
    $ bqueues    # available queues
    $ lsload     # check load values for all hosts

In class exercises
==================
::
 1. test
 2. test 2
