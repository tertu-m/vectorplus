VectorPLUS
==========
VectorPLUS is a bioinformatics aide written by prototypical college students Henry Ward and Jack Walstrom.

Currently, it can do one thing and one thing only: find a given gene sequence within a larger sequence and
outputs its regulatory region, given the size of said region. Reads and writes using fasta files.

USE: Put your larger DNA sequence, likely a scaffold, and your search sequence in the same folder as
VectorPLUS.py. Make sure they are fasta files. Then, run VectorPLUS. 

UPCOMING FEATURES: Fuzzy searching (finally)

KNOWN ISSUES: Not sufficiently tested yet for cases where multiple equal size matches are found. 
Has a minimum search size of 5. 
