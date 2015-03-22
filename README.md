
#Oases
_De novo_ transcriptome assembler based on the Velvet genome assembler core

##Requirements
Oases should function on any standard 64 bit Unix environment with a C compiler.
A good amount of physical memory (12GB to start with, more is no luxury) is recommended.

##Installing Oases

    > git clone --recursive https://github.com/dzerbino/oases
    > cd oases
    > make
    > ./oases --version
    1.2.08

##Compilation options

There are various compile time options that can be used:

| Option | Default | Description 
| ------ | ------- | ------------
| CATEGORIES | 2 | Maxium number of different DNA libraries
| MAXKMERLENGTH | 64 | Maximum k-mer value supported
| OPENMP | 0 | Enable OpenMP multithreading support
| BIGASSEMBLY | 0 | FIXME
| VBIGASSEMBLY | 0 | FIXME
| LONGSEQUENCES | 0 | FIXME
| SINGLE_COV_CAT | 0 | FIXME
| BUNDLEDZLIB | 0 | Use the bundled zlib 1.2.3 instead of system one

You can apply them as in this example which changes three options:

    make OPENMP=1 CATEGORIES=4 MAXKMERLENGTH=192


#Running Oases

If you are not familiar with using Velvet, read Velvet’s manual first,
as many references below will point to that document.

##For impatient people

    > velveth directory 21,23 data/test_reads.fa
    > velvetg directory_21 -read_trkg yes
    > oases directory_21
    > ls directory_21

    > velvetg directory_23 -read_trkg yes
    > oases directory_23

    > velveth mergedAssembly 23 -long directory*/transcripts.fa
    > velvetg mergedAssembly -read_trkg yes -conserveLong yes
    > oases mergedAssembly -merge yes

Or use the python script ```oases_pipeline.py```:

    > python scripts/oases_pipeline.py -m 21 -M 23 -d 'data/test_reads.fa'

##For patient people

To run Oases it is recommended to run an array of single-_k_ assemblies,
for different values of _k_ (i.e. the hash length). These assemblies are
then merged into a final assembly.

###Single-_k_ assemblies

Each single-_k_ assembly consists of a simple Velvet run, followed by an
Oases run. As in Velvet, the hash length _k_ is an odd integer. Velveth
is run normally and velvetg is run with only one option:

    > velveth directory_k k reads.fa
    > velvetg directory_k -read_trkg yes

If you have strand specific reads then you have to type:

    > velveth directory_k k -strand_specific reads.fa

Oases is then run on that directory:

    > oases directory_$k

###Oases options

### Oases help

In case of doubt, you can always print a help message by typing
(interchangeably):

    > oases
    > oases --help

### Insert length

If using paired-end reads, Oases needs to know the expected insert
length (it cannot be estimated reliably automatically). The syntax to
provide insert lengths is identical to that used in Velvet:

    oases directory -ins_length 500

### Coverage cutoff

Just like Velvet, low coverage contigs are removed from the assembly.
Because genes span all the spectrum of expression levels and therefore
coverage depths, setting the coverage cutoff does not depend on an
expected or median coverage. Instead, the coverage cutoff is set to
remove all the contigs whose coverage is so low that *de novo* assembly
would be impossible anyways (by default it is set at $3x$). If you wish
to set it to another value, simply use the Velvet-like option:

    oases directory -cov_cutoff 5

### Edge fraction cutoff

To allow more efficient error removal at high coverage, Oases integrates
a dynamic correction method: if at a given node, an edge’s coverage
represents less than a given percentage of the sum of coverages of edges
coming out of that node, it is removed.

By default this percentage is 10% but you can set it manually:

    oases directory -edgeFractionCutoff 0.01

### Scaffold filtering

By default the distance estimate between two long contigs is discarded
as unreliable if it is supported by less than a given number (by default
4) or read-pairs or bridging reads. If you wish to change this cutoff:

    oases directory -min_pair_count 5

### Filtering the output

By setting the minimum transfrag length (by default 100bp), the user can
control what is being output in the results files:

    > oases directory -min_trans_lgth 200

### Assembly merging

After running the previous process for different values of _k_ it is
necessary to merge the results of all these assemblies (contained in
```transcripts.fa``` files) into a single, non redundant assembly. To
realize this merge, it is necessary to choose a value of $K$.
Experiments have shown that $K=27$ works nicely for most organisms.
Higher values for $K$ should be used if problems with missassemblies are
observed. 

Assume we want to create a merged assembly in a folder called
MergedAssembly.

    > velveth MergedAssembly/ K -long directory*/transcripts.fa
    > velvetg MergedAssembly/ -read_trkg yes -conserveLong yes
    > oases MergedAssembly/ -merge yes

 that the transcripts.fa files need to be given as *long*.\

### Using the supplied python script

With version 0.2 of Oases comes a new python script that runs the
individual single-_k_ assemblies and the new merge module. We highly
recommend to use the script for the analyses, but please read the
previous subsections [insert]-[filtering] for the Oases params that you
need to supply to the script. However, as the script also runs velvet
please consult the velvet manual for more insight. First notice the
options for the script below

    python scripts/oases_pipeline.py
    Options:
      -h, --help            show this help message and exit
      -d FILES, --data=FILES
                            Velveth file descriptors
      -p OPTIONS, --options=OPTIONS
                            Oases options that are passed to the command line,
                            e.g., -cov_cutoff 5 -ins_length 300
      -m KMIN, --kmin=KMIN  Minimum k
      -M KMAX, --kmax=KMAX  Maximum k
      -s KSTEP, --kstep=KSTEP
                            Steps in k
      -g KMERGE, --merge=KMERGE
                            Merge k
      -o NAME, --output=NAME
                            Output directory prefix
      -r, --mergeOnly       Only do the merge
      -c, --clean           Clean temp files

 that the script checks if Oases version at least 0.2.01 and Velvet at
least 1.1.7 are installed on your path. Otherwise it will not work.

We will illustrate the usage of the script with the two most common
scenarios single-end and paired-end reads.\
First assume we  want to compute a merged assembly from 21 to 35 with
Oases on a single-end data set that is strand specific and retain
transcripts of minimum length 100. We want to save the assemblies in
folders that start with singleEnd :

    python scripts/oases_pipeline.py -m 21 -M 35 -o singleEnd 
    -d ' -strand_specific yes data/test_reads.fa '
    -p ' -min_trans_lgth 100 ' 

The script creates a folder named singleEnd\__k_ for each single-_k_
assembly and one folder named singleEndMerged that contains the merged
transcripts for all single-_k_ assemblies in the *transcripts.fa* file.

Now assume we have paired-end reads with insert length 300 and we want
to compute a merged assembly from 17 to 40. We want to save the
assemblies in the folders that start with pairedEnd

    python scripts/oases_pipeline.py -m 17 -M 40 -o pairedEnd 
    -d ' -shortPaired  data/test_reads.fa '
    -p ' -ins_length 300 '

Now say we found that using k=17 was a bit too low and we want to look
at a different merged assembly range. Instead of redoing all the time
consuming assemblies we can re-run the merged module on a different
range from 21 to 40 like this:

    python scripts/oases_pipeline.py -m 21 -M 40 -r -o pairedEnd  

  that it is essential when you re-run the merged assembly for the
script to function properly to give the exact same output prefix with
the $-o$ parameter.

## Output files

Oases produces two output files. The predicted transcript sequences can
be found in *transcripts.fa* and the file *contig-ordering.txt* explains
the composition of these transcripts, see below.

1.  new_directory/transcripts.fa – A FASTA file containing the
    transcripts imputed directly from trivial clusters of contigs (loci
    with less than two transcripts and Confidence Values = 1) and the
    highly expressed transcripts imputed by dynamic programming (loci
    with more than 2 transcripts and Confidence Values $<$1).

2.  new_directory/contig-ordering.txt – A hybrid file which describes
    the contigs contained in each locus in FASTA format, interlaced with
    one line summaries of the transcripts generated by dynamic
    programming. Each line is a string of atoms defined as:
    \$contig\_id:\$cumulative\_length-(\$distance\_to\_next\_contig)-$>$next
    item

    Here the cumulative length is the total length of the transcript
    assembly from its 5’ end to the 3’ end of that contig. This allows
    you to locate the contig sequence within the transcript sequence.

## FAQ and Practical considerations

### I am running out of memory. What should I do? 

Depending on the complexity of the dataset it may happen that Velvet and
Oases start to use a large amount of memory. In general, the assembly
graphs for transcriptome data are smaller as compared to genome data.
However, because the coverage is not even roughly uniform it depends on
the sample being sequenced how much memory is used. Below we give
examples of pre-processing steps that lead to a decrease in memory
consumption, sorted in order of importance.

####Avoid overlapping reads

In our experience, if insert lengths are so short that paired reads
overlap on their 3’ ends, i.e., if insert length $<$ 2 $\cdot$ read
length, this prevents the early filtering of low coverage reads. This
means that even small datasets may require huge amounts of memory. A way
to avoid the problem is to join paired-end reads that overlap for
example using the software FLASH *(FLASH: fast length adjustment of
short reads to improve genome assemblies, Magoc T, Salzberg SL.
Bioinformatics 2011)*. Otherwise it is possible to clip reads so that
their overlap is strictly shorter than the hash length, that might lead
to a huge loss of data however.\

#### Avoid bad quality reads

Reads with many errors create new nodes in the de Bruijn graph, all of
which are later discarded by the error removal steps. Removing the
errors in the reads *a priori* is a good way to decrease the memory
consumption. A very common approach is to remove bad quality bases from
the ends of read sequences or remove entire read sequences, numerous
packages exist to do that but differ in functionality and definition of
bad quality (in no particular order):

1.  FASTX-toolkit - <http://hannonlab.cshl.edu/fastx_toolkit/> (diverse
    functionality, command line, Galaxy support)

2.  ea-utils - <http://code.google.com/p/ea-utils/wiki/FastqMcf>
    (diverse functionality, command line)

3.  ShortRead (R-package) - see here for long version
    <http://manuals.bioinformatics.ucr.edu/home/ht-seq> or short version
    <http://jermdemo.blogspot.com/2010/03/soft-trimming-in-r-using-shortread-and.html>

4.  Condetri - <http://code.google.com/p/condetri/> (only read trimmer
    good documentation, command line)

5.  bwa - <http://bio-bwa.sourceforge.net/bwa.shtml> (command line *aln*
    option)

Note that when you use paired-end sequences and reads are discarded from
the file you need to make sure that the fasta/fastq file that you give
to velveth contains the paired-end reads in correct pairing. Also if the
read clipper retains single end reads from a pair, it is advisable to
give them as a second dataset of type *-short* to avoid data loss.

#### Trim sequencing adapters

Depending on your sequencing experiment, some of the reads may have
adapters that were used, e.g, during library preparation. These adapters
should be trimmed. You can use the Fastx-toolkit or ea-utils to do that,
among many others, see above.

#### Remove duplicate reads

Especially for paired-end datasets it happens that the exact same read
sequence is found a large number of times. Again you can use some of the
tools mentioned above to remove these duplicates.

### Oases does not give an expression level. How can I get that? 

In general you need to align reads to the transcripts with your favorite
aligner, one that allows for multi mappings !! Although few results
exist about the difficulties and problems of quantitation of *de novo*
transcripts’ expression levels, here are a few pointers to specialized
tools:

1.  eXpress - <http://bio.math.berkeley.edu/eXpress/> (command line,
    Windows support)

2.  RSEM - <http://deweylab.biostat.wisc.edu/rsem/> (command line)

### I want to predict coding regions. How can I do that? 

All the transcripts reported in a locus should have the same
directionality, but that is only true if there is no contamination.
Standard approaches are to use blastp and tools alike. Alternatively use
gene prediction tools like GlimmerHMM, SNAP or AUGUSTUS.

## Mailing list:

For questions/requests/etc. you can subscribe to the users’ mailing list: 
```oases-users@ebi.ac.uk```. 
You can sign up at the
[oases-users listserver](http://listserver.ebi.ac.uk/mailman/listinfo/oases-users).

##Paper
Schulz MH, Zerbino DR, Vingron M, Birney E. 
_Oases: robust de novo RNA-seq assembly across the dynamic range of expression levels._
__Bioinformatics__ 2012 Apr 15;28(8):1086-92. 
[PubMed](http://www.ncbi.nlm.nih.gov/pubmed/22368243)

##Authors
* Daniel Zerbino <dzerbino@soe.ucsc.edu>
* Marcel Schulz <marcel.schulz@molgen.mpg.de>

