*************************
Programming Tips & Tricks
*************************

Overview
========
Several things influence how effectively you learn to program, and how
well you program once you have mastered the basic ideas. A major issue is
your efficiency in actaully *using* a computer, and not programming *per
se*.

For example, the longer you spend searching for a particular key to type,
or surfing around with your mouse, the less time you spend writing,
running and debugging programs. Here are several pointers to help you
become more efficient at using computers, independent of learning
programming languages.

Learn to type
=============
Hunting and pecking is inefficient, and prevents you from spending your
valuable time efficiently. If you're looking at your keyboard, you're
*not* looking at the screen, reading and debugging code. Once you have the
typing basics down, you should be typing 40-50 words per minute, without
ever looking at the keyboard.

There are several good, modern tools [#]_ to help you master touch-typing.

.. [#] 
    - Touch Typing Tutorial : http://www.typingweb.com/
    - Typing IO : Language specific typing practice http://typing.io/

Learn to type funny characters
------------------------------
In programming you use a variety of characters that you don't always use
typing other kinds of documents. Learn the locations of the following by
heart:

+ Number sign (for commenting): #
+ Dollar sign (for variables): $
+ Underscore (for variable naming): _
+ Parentheses: ()
+ Brackets: []
+ Curly brackets: {}
+ Tilde (i.e. the squiggle; for going $HOME): ~
+ Math symbols ("+", "-", "*", "/", "=")

Learn hot keys for window management
====================================
**The mouse is your enemy.** Yes, it revolutionized the computer–human
interaction. But the more time you spend using your mouse, the less time
you spend with your hands on the keyboard and doing useful things.

You can do most things with your keyboard. There are several hot keys you
should learn that will maximize your productivity on the computer by
minimizing your use of the mouse:

+ **<Alt>-<Tab>** : Flip through windows quickly and effortlessly
  without ever touching your mouse.
+ **<Ctrl>-<Page Up/Down>** : Switch between Terminal windows on Linux

.. tip::

    **Launch your most-used apps automatically during login.**

    For example, automatically launch four terminal windows and a browser
    window, without having to click.

Learn to use a terminal text editor
===================================
There are two types of nerds in this world: 

- ``vim`` users
- ``emacs`` users
    
I'm a ``vim``-user. I can't even log out of an ``emacs`` session.

Learning a terminal text editor like ``vim`` increases productivity
substantially, because it allows you to:

+ run the editor within an existing terminal, without opening a new
  window
+ work on multiple documents simultaneously
+ syntax highlight your code
+ manipulate blocks of text quickly and efficiently

You can run ``vim`` from the terminal prompt:

.. code-block:: bash

    $ vim filename.txt

To quit a vim session, you need to:

#. enter `command mode` with the colon key
#. write the file
#. quit the program

This can be accomplished by typing::

    :wq <enter>

In your copious spare time, and after you have mastered the basics of
shell, Python and R programming, you should take a tutorial [#]_ on using
a terminal text editor.

.. [#] OpenVim http://www.openvim.com/ 

.. raw:: pdf

    PageBreak
