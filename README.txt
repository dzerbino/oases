README.TXT

OASES SOURCE
Feb 1, 2010
Daniel Zerbino (dzerbino@soe.ucsc.edu)
Marcel Schulz (marcel.schulz@molgen.mpg.de)

> SUMMARY
        * A/ REQUIREMENTS
        * B/ COMPILING INSTRUCTIONS
	* C/ RUNNING

----------------------------------------------------------------------------------
A/ REQUIREMENTS

        Oases should function on any standard 64bit Linux environment with
gcc. A good amount of physical memory (12GB to start with, more is no luxury)
is recommended.
	Before trying to compile Oases, you must install the Velvet package:
www.ebi.ac.uk/~zerbino/velvet/ . Keep note of the directory in which you install
Velvet.

----------------------------------------------------------------------------------
B/ COMPILING INSTRUCTIONS

Normally, with a GNU environment, just type:

> make 'VELVET_DIR=/path/to/velvet'

Note that you need to communicate all the Velvet compilation settings during
the Oases compilation. Therefore, if you want to make a debugging colorspace version of
Oases with a maximum kmer length of 63 and 5 short-read libraries, the
commandline becomes:

> make colordebug 'VELVET_DIR=/path/to/velvet' 'MAXKMERLENGTH=63'\
'CATEGORIES=5'

----------------------------------------------------------------------------------
C/ RUNNING

The compilation should have produced a file names OasesManual.pdf with all you 
could want to know about running Oases. 
