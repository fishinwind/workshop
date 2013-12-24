Class 3: The command-line++
===========================

Goals
-----

1. understand how to apply some common linux utilities to files
2. understand pipes (|)


Compressed Files
----------------

Most common way to compress single files is `gzip`

    $ ls /opt/bio-workshop/data/\*.gz

We can (g)unzip that file as:

    $ gunzip /opt/bio-workshop/data/some.gz

And re-zip is as:

    $ gzip /opt/bio-workshop/data/some

But if we just want to stream the uncompressed data without changing the file:

    $ zless /opt/bio-workshop/data/some.gz

Pipes
-----

We probably want to do something with the file as we uncompress it:

    $ zless /opt/bio-workshop/data/some.gz | head

We already know the head command prints the first -n lines.

Try piping the output to some other commands (tail|echo|cowsay)


