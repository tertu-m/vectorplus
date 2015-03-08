VectorPLUS (0.25)
==========
VectorPLUS is a bioinformatics aide written by prototypical college students Henry Ward and Jack Walstrom.

VectorPlus returns regions upstream of certain positions in a DNA/RNA sequence.
These positions can be determined in one of two ways:
1) Position input - enter the desired positions to find regions upstream of
2) Sequence search - enter a smaller sequence to search for in a larger sequence

Reads and writes using FASTA files. Works best if the FASTA file you search in has only one sequence, such
as a scaffold.

NOTE: Make sure all of your sequences to be searched in/with are in the same folder as VectorPlus. Sequence search is a rather dumb algorith, and despite certain capabilities (essentially replaces search capabilities in VectorNTI) will soon be replaced with local BLAST. 

UPCOMING FEATURES: BLAST searching using BioPython. GUI to visualize sequences. 

