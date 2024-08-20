==================
counting_sequences
==================


This contains a simple sequence counter intended to count the occurences of 
a set of predefined sequences in a fastq file.



Description
===========

This contains a simple sequence counter intended to count the occurences of 
a set of predefined sequences in a fastq file. First, all occuring read sequences
are counted. Then, these counts are mapped to a set of predefined sequences.

It also contains a wrapper for NGmerge to merge overlapping paired-end
sequenced reads.

Reads and sequences can also be trimmed to exclude flanking regions that
are not of interest, e.g. adapters and primers.



Note
====

This project is intended for use with pypiegraph2 and mbf packages available
on pypi and github for job generation and scheduling.
