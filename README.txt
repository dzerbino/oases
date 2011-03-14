README.TXT

OASES SOURCE
Feb 1, 2010
Daniel Zerbino (dzerbino@soe.ucsc.edu)
Marcel Schulz (marcel.schulz@molgen.mpg.de)

> SUMMARY
        * A/ REQUIREMENTS
        * B/ COMPILING INSTRUCTIONS
	* C/ RUNNING
	* D/ OUTPUT FILES
	* E/ OPTIONS

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

You must first process the reads using Velvet: 
* you must choose a hash length at this stage (cf. the Velvet manual),
* DO NOT set a coverage cutoff, you should set that when running oases, 
* DO NOT set an expected coverage,
* remember to turn on the -read_trkg option when running velvetg. 

As an example:

> velveth new_directory 21 -shortPaired data/test_reads.fa
> velvetg new_directory -read_trkg yes

You can now run Oases on the Velvet working directory which has just been created.
Provide all the information about insert lengths and their standard deviation as 
possible (identical to Velvet):

> oases new_directory -ins_length 200

----------------------------------------------------------------------------------
D/ OUTPUT FILES

Oases produces a number of output files, which correspond to the different algorithms
being run succesively on the data. In the above example, you would find:

new_directory/transcripts.fa
	A FASTA file containing the transcripts imputed directly from trivial
	clusters of contigs (loci with less than two transcripts and Confidence Values = 1)
	and the highly expressed transcripts imputed by dynamic
	programming (loci with more than 2 transcripts and Confidence Values <1).

new_directory/splicing_events.txt
	A hybrid file which describes the contigs contained in each locus in FASTA
	format, interlaced with one line descriptions of splicing events using the
	AStalavista nomenclature*.

new_directory/contig-ordering.txt
	A hybrid file which describes the contigs contained in each locus in FASTA
	format, interlaced with one line summaries of the transcripts listed
	in transcripts.fa . Each line is a string of atoms defined as:
	$contig_id:$cumulative_length-($distance_to_next_contig)->

	Here the cumulative length is the total length of the transcript assembly from
	its 5' end to the 3' end of that contig. This allows you to locate the contig
	sequence within the transcript sequence.

	A file describing the transcripts imputed directly from trivial clusters.
	Its format is identical to the file described previously.


* Sammeth, Michael, Foissac, Sylvain  Guigó, Roderic, 'A General Definition and 
Nomenclature for Alternative Splicing Events', PLoS Comput Biol , vol. 4, no. 8, 
e100014y+ (2008). 

----------------------------------------------------------------------------------
E/ OPTIONS

The behavior of Oases can be modified using the following options:

-min_trans_length
	simple threshold on output transcript length
-cov_cutoff 
	minimum number of times a k-mer has to be observed to be used in the 
	assembly (just like in Velvet) [default=3]
-min_pair_cov
	minimum number of times two contigs must be connected by reads or read pairs
	to be clustered together [default=4]
-paired_cutoff 
	minimum ratio between the numbers of observed and expected connecting
	read pairs between two contigs [default=0.1]
-scaffolding
	allows you to prevent the creation of gapped transcripts

E.g.:

> oases new_directory -ins_length 200 -cov_cutoff 3 -min_pair_count 4

Running:
> oases --help 
will produce a short help message

