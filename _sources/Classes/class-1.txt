.. include:: /_static/substitutions.txt

*****************************************************
             Class 1 : Class Introduction
*****************************************************

:Class date: |c1-date| 

Goals
=====
#. Class overview
#. Get the VM running
#. Overall goals for the Classes

Class Overview
==============
Each class is 2 hours. We intend to spend the first 60 min going
through exercises that demonstrate how specific tools are useful in
bioinformatics. During the remaining hour, we expect you to work through
exercises, asking for help when you get stuck. 

We will record the first 60 minutes using Panopto Screen Capture, and
these recordings will be available in Canvas. We have found that simply
watching someone work in a terminal (move around, open up text editors,
write and execute simple programs) can be a very effective way to get
started with programming. Hopefully these movies will serve the same
purpose.

Each week, we will have 1 take home quiz, due the following Tuesday at 5
PM. 

Linux installations
===================
The PCs in the library have Virtual Box installed with a minimal Linux
installation. If you have your own PC laptop, you can install Virtual Box
and any standard Linux distribution (Ubuntu or Mint). If you have a Mac
laptop, you can do the same, or just use the native terminal.

In any case, we will create logins on our compute cluster (amc-tesla) and
all your work will be done through that.

Terminal and text editors
=========================
When you open a Terminal, you also launch a shell process, typically a
bash process. At the prompt, you can type things that bash understands,
and it will do them. The shell has its own language, which you will learn
over the course of the class. It also runs executable files that it can
find on the PATH (i.e. the set of directories that contain exectuables).

You can find what is on your PATH by typing:

.. code-block:: bash

   $ echo $PATH

The PATH is one of several environment variables that are created when you
login. You can see all of these with:

.. code-block:: bash

   $ env

.. nextslide::
    :increment:

One important program in the PATH is `vim`. You will use this program to
keep notes and write small programs. 

You can run ``vim`` from the terminal prompt:

.. code-block:: bash

    $ vim filename.txt

You should notice that the prompt will disappear and you will be in a
`vim` session.

Now press the `i` key to enter `insert` mode, and start typing. Press
`ESC` to exit `insert` mode.

To quit a vim session, you need to:

#. enter `command mode` with the colon key
#. write the file
#. quit the program

This can be accomplished by typing::

    :wq <enter>

Practice using `vim` with this tutorial [#]_.

.. [#] OpenVim http://www.openvim.com/ 

Cluster access
==============
We have set up accounts for the class on our departmental cluster. We will
set up your accounts at the end of class and reset your passwords:

.. code-block:: bash

    # the -X flag starts an X11 connection 
    $ ssh -X username@amc-tesla.ucdenver.pvt

    ...

    # once you are logged in, text your X11 connection with
    $ xeyes

Running jobs on the cluster
===========================
First you will grab a single CPU from the queueing system so that you can
work without affecting the head node. We use ``qlogin`` for this:

.. code-block:: bash

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

.. nextslide::
    :increment: 


Shell and Python Programming
============================
It is important that you learn a few new computer languages. Others have
developed very good guides to teach you these languages, and we are going
to use those in the class. We expect you to begin taking these classes
immediately.

You will spend a lot of time going through these online classes, both in
scheduled class time, and outside of class time. Instead of focusing on
teaching you these languages, we will focus on helping you get through all
of the frustating problems that come up when you're learning the languages.

We will spend the first ~2 weeks learning shell [#]_ and all the things
you have access to within the shell.

.. [#] The Command Line Crash Course
        http://cli.learncodethehardway.org/book/

After learning the shell, we will begin learning R and several packages
within R.

Finally, we will begin learning Python [#]_. The Python language allows
you to do more sophisticated things that would be possible in shell or R, but
would be considerably more clunky.

.. [#] Learn Python the Hard Way
        http://learnpythonthehardway.org/book/

First Quiz : Reading
====================
Computational biology projects inevitably accrue a lot of files. For the
first quiz, you'll need to read a paper [#]_ and be able to put a set of
files in the correct places. We require that you use this scheme for
all of your projects in and out of the class.

.. [#] A Quick Guide to Organizing Computational Biology Projects (2009)
        PLoS Comput. Biol. William S. Noble
        http://dx.plos.org/10.1371/journal.pcbi.1000424

.. raw:: pdf

    PageBreak

