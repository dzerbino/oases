/*
    Copyright 2009,2010 Daniel Zerbino (dzerbino@soe.ucsc.edu)

    This file is part of Oases.

    Oases is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    Oases is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Oases. If not, see <http://www.gnu.org/licenses/>.
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// Compilation
#include "globals.h"

// Utilities
#include "graphStats.h"
#include "utility.h"

// Datastructures
#include "kmer.h"
#include "readSet.h"
#include "tightString.h"
#include "graph.h"

// Graph operations
#include "graph.h"
#include "graphReConstruction.h"
#include "concatenatedGraph.h"
#include "correctedGraph.h"

// Oases specific stuff 
#include "transcript.h"

static int OASES_VERSION_NUMBER = 0;
static int OASES_RELEASE_NUMBER = 1;
static int OASES_UPDATE_NUMBER = 12;

static void printUsage()
{
	puts("Usage:");
	puts("./oases directory [options]");
	puts("");
	puts("\tdirectory\t\t\t: working directory name");
	puts("");
	puts("Standard options:");
	puts("\t-ins_length2 <integer>\t\t: expected distance between two paired-end reads in the second short-read dataset (default: no read pairing)");
	puts("\t-ins_length_long <integer>\t: expected distance between two long paired-end reads (default: no read pairing)");
	puts("\t-ins_length*_sd <integer>\t: est. standard deviation of respective dataset (default: 10% of corresponding length)");
	puts("\t\t[replace '*' by nothing, '2' or '_long' as necessary]");
	puts("\t-unused_reads <yes|no>\t\t: export unused reads in UnusedReads.fa file (default: no)");
	puts("\t-amos_file <yes|no>\t\t: export assembly to AMOS file (default: no export)");
	puts("\t-alignments <yes|no>\t\t: export a summary of contig alignment to the reference sequences (default: no)");
	puts("\t--help\t\t\t\t: this help message");
	puts("Advanced options:");
	puts("\t-cov_cutoff <floating-point>\t: removal of low coverage nodes AFTER tour bus or allow the system to infer it (default: 3)");
	puts("\t-min_pair_count <integer>\t: minimum number of paired end connections to justify the scaffolding of two long contigs (default: 4)");
	puts("\t-min_trans_lgth <integer>\t: Minimum length of output transcripts (default: hash-length)");
	puts("\t-paired_cutoff <floating-point>\t: minimum ratio allowed between the numbers of observed and estimated connecting read pairs");
	puts("\t\tMust be part of the open interval ]0,1[ (default: 0.1)");
	puts("\t-scaffolding <yes|no>\t\t:Allow gaps in transcripts (default: yes)");
	puts("\t-degree_cutoff <integer>\t: Maximum allowed degree on either end of a contigg to consider it 'unique' (default: 3)");
	puts("");
	puts("Output:");
	puts("\tdirectory/transcripts.fa");
	puts("\tdirectory/splicing_events.txt");
	puts("\tdirectory/contig-ordering.txt");
}

int main(int argc, char **argv)
{
	Graph *graph;
	char *directory, *graphFilename, *seqFilename, *transcriptFilename,
	    *eventFilename;
	Coordinate insertLength[CATEGORIES];
	Coordinate insertLengthLong = -1;
	Coordinate std_dev[CATEGORIES];
	Coordinate std_dev_long = -1;
	FILE *file;
	int arg_index, arg_int;
	char *arg;
	IDnum *sequenceLengths = NULL;
	Category cat;
	long long longlong_var;
	short int short_var;
	Locus *loci;
	IDnum locusCount;
	ReadSet *reads;
	double coverageCutoff = 3;
	boolean *dubious = NULL;
	Coordinate minTransLength = 0;
	double pairedThreshold = 0.1;
	boolean unusedReads = false;
	boolean exportAssembly = false;
	boolean scaffolding = true;
	boolean exportAlignments = false;

	setProgramName("oases");

	for (cat = 0; cat < CATEGORIES; cat++) {
		insertLength[cat] = -1;
		std_dev[cat] = -1;
	}

	// Error message
	if (argc == 1) {
		puts("oases - De novo transcriptome assembler for the Velvet package");
		printf("Version %i.%i.%2.2i\n", OASES_VERSION_NUMBER,
		       OASES_RELEASE_NUMBER, OASES_UPDATE_NUMBER);
		puts("\nCopyright 2009,2010 Daniel Zerbino (dzerbino@soe.ucsc.edu)");
		puts("This is free software; see the source for copying conditions.  There is NO");
		puts("warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
		puts("Compilation settings:");
		printf("CATEGORIES = %i\n", CATEGORIES);
		printf("MAXKMERLENGTH = %i\n", MAXKMERLENGTH);
		puts("");
		printUsage();
		return 1;
	}

	if (strcmp(argv[1], "--help") == 0) {
		printUsage();
		return 0;
	}
	// Memory allocation 
	directory = argv[1];
	graphFilename = mallocOrExit(strlen(directory) + 100, char);
	seqFilename = mallocOrExit(strlen(directory) + 100, char);
	transcriptFilename = mallocOrExit(strlen(directory) + 100, char);
	eventFilename = mallocOrExit(strlen(directory) + 100, char);

	// Argument parsing
	for (arg_index = 2; arg_index < argc; arg_index++) {
		arg = argv[arg_index++];
		if (arg_index >= argc) {
			puts("Unusual number of arguments!");
			printUsage();
			exit(1);
		}
		if (strcmp(arg, "-cov_cutoff") == 0) {
			sscanf(argv[arg_index], "%lf", &coverageCutoff);
		} else if (strcmp(arg, "-min_trans_lgth") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			minTransLength = (Coordinate) longlong_var;
		} else if (strcmp(arg, "-degree_cutoff") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setDegreeCutoff(arg_int);
		} else if (strcmp(arg, "-paired_cutoff") == 0) {
			sscanf(argv[arg_index], "%lf", &pairedThreshold);
			setPairedThreshold(pairedThreshold);
		} else if (strcmp(arg, "-alignments") == 0) {
			exportAlignments =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-amos_file") == 0) {
			exportAssembly =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-unused_reads") == 0) {
			unusedReads =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-scaffolding") == 0) {
			scaffolding =
			    (strcmp(argv[arg_index], "yes") == 0);
		} else if (strcmp(arg, "-min_pair_count") == 0) {
			sscanf(argv[arg_index], "%i", &arg_int);
			setUnreliableConnectionCutoff_oases(arg_int);
		} else if (strcmp(arg, "--help") == 0) {
			printUsage();
			return 0;
		} else if (strcmp(arg, "-ins_length") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLength[0] = (Coordinate) longlong_var;
			if (insertLength[0] < 0) {
				printf("Invalid insert length: %lli\n",
				       (long long) insertLength[0]);
				exit(1);
			}
		} else if (strcmp(arg, "-ins_length_sd") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev[0] = (Coordinate) longlong_var;
			if (std_dev[0] < 0) {
				printf("Invalid std deviation: %lli\n",
				       (long long) std_dev[0]);
				exit(1);
			}
		} else if (strcmp(arg, "-ins_length_long") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLengthLong = (Coordinate) longlong_var;
		} else if (strcmp(arg, "-ins_length_long_sd") == 0) {
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev_long = (Coordinate) longlong_var;
		} else if (strncmp(arg, "-ins_length", 11) == 0
			   && strchr(arg, 'd') == NULL) {
			sscanf(arg, "-ins_length%hi", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				printf("Unknown option: %s\n", arg);
				exit(1);
			}
			sscanf(argv[arg_index], "%lli", &longlong_var);
			insertLength[cat - 1] = (Coordinate) longlong_var;
			if (insertLength[cat - 1] < 0) {
				printf("Invalid insert length: %lli\n",
				       (long long) insertLength[cat - 1]);
				exit(1);
			}
		} else if (strncmp(arg, "-ins_length", 11) == 0) {
			sscanf(arg, "-ins_length%hi_sd", &short_var);
			cat = (Category) short_var;
			if (cat < 1 || cat > CATEGORIES) {
				printf("Unknown option: %s\n", arg);
				exit(1);
			}
			sscanf(argv[arg_index], "%lli", &longlong_var);
			std_dev[cat - 1] = (Coordinate) longlong_var;
			if (std_dev[cat - 1] < 0) {
				printf("Invalid std deviation: %lli\n",
				       (long long) std_dev[cat - 1]);
				exit(1);
			}
		} else {
			printf("Unknown option: %s;\n", arg);
			printUsage();
			return 1;
		}
	}

	// Bookkeeping
	logInstructions(argc, argv, directory);

	strcpy(seqFilename, directory);
	strcat(seqFilename, "/Sequences");

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/Graph2");

	// Graph uploading or creation
	if ((file = fopen(graphFilename, "r")) != NULL) {
		fclose(file);
		strcpy(graphFilename, directory);
		strcat(graphFilename, "/Graph2");
		graph = importGraph(graphFilename);
		reads =
		    importReadSet(seqFilename);
		convertSequences(reads);
	} else {
		puts("No Graph2 file to work with!");
		puts("Please re-run Velvetg with the -read_trkg option on.");
		return 1;
	}

	sequenceLengths =
	    getSequenceLengths(reads, getWordLength(graph));

	correctGraph(graph, sequenceLengths, reads->categories);

	dubious =
	    removeLowCoverageNodesAndDenounceDubiousReads(graph,
							  coverageCutoff,
							  reads);
	clipTipsHard(graph);

	// Set insert lengths and their standard deviations
	createReadPairingArray(reads);
	for (cat = 0; cat < CATEGORIES; cat++) {
		if (insertLength[cat] > -1 && std_dev[cat] < 0)
			std_dev[cat] = insertLength[cat] / 10;
		setInsertLengths(graph, cat,
				 insertLength[cat], std_dev[cat]);
		if (insertLength[cat] > -1) {
			pairUpReads(reads, 2 * cat + 1);
		}

		detachDubiousReads(reads, dubious);
	}

	if (insertLengthLong > -1 && std_dev_long < 0)
		std_dev_long = insertLengthLong / 10;
	setInsertLengths(graph, CATEGORIES, insertLengthLong,
			 std_dev_long);
	if (insertLengthLong > -1)
		pairUpReads(reads, 2 * CATEGORIES + 1);

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/stats.txt");
	displayGeneralStatistics(graph, graphFilename, reads);

	strcpy(graphFilename, directory);
	strcat(graphFilename, "/LastGraph");
	exportGraph(graphFilename, graph, reads->tSequences);

	loci =
	    extractGraphLoci(graph, reads, dubious, sequenceLengths, &locusCount, scaffolding);

	computeTranscripts(loci, locusCount);
	strcpy(transcriptFilename, directory);
	strcat(transcriptFilename, "/transcripts.fa");
	exportTranscripts(loci, locusCount, transcriptFilename, minTransLength);
	strcpy(transcriptFilename, directory);
	strcat(transcriptFilename,
	       "/contig-ordering.txt");
	exportContigOrders(loci, locusCount, transcriptFilename, minTransLength);
	printf("Finished heuristic approach, used %li/%li reads\n", (long) usedTranscriptReads(graph, minTransLength, loci, locusCount), (long) sequenceCount(graph));

	if (unusedReads) 
		exportUnusedTranscriptReads(graph, loci, locusCount, reads, minTransLength, directory);

	if (exportAssembly) 
		exportAMOSTranscripts(graph, loci, locusCount, reads, minTransLength, directory);

	if (exportAlignments)
		exportTranscriptMappings(loci, locusCount, graph, reads, minTransLength, directory);

	cleanTranscriptMemory(loci, locusCount);

	removeIndirectConnections();

#ifndef COLOR
	computeASEvents(loci, locusCount);
	strcpy(eventFilename, directory);
	strcat(eventFilename, "/splicing_events.txt");
	exportASEvents(loci, locusCount, eventFilename);
#endif

	cleanLocusMemory(loci, locusCount);
	destroyGraph(graph);
	free(graphFilename);
	free(seqFilename);
	free(transcriptFilename);
	free(eventFilename);
	free(dubious);
	destroyReadSet(reads);
	return 0;
}
