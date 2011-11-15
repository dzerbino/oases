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
#ifndef _TRANSCRIPT_H_
#define _TRANSCRIPT_H_

typedef struct locus_st Locus;
typedef struct transcript_st Transcript;
typedef struct event_st Event;

struct locus_st {
	IDnum contigCount;
	IDnum longContigCount;
	Node **contigs;
	Transcript *transcript;
	Event *event;
};

struct transcript_st {
	IDnum contigCount;
	Node **contigs;
	Coordinate *distances;
	double confidence;
	Transcript *next;
};

void setUnreliableConnectionCutoff_oases(int val);

Locus *extractGraphLoci(Graph * graph, ReadSet * argReads,
			boolean * dubious, ShortLength * lengths,
			IDnum * locusCount, boolean scaffolding);
void computeTranscripts(Locus * loci, IDnum locusCount);
void computeASEvents(Locus * loci, IDnum locusCount);

void exportTranscripts(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength, Graph * graph);
void exportASEvents(Locus * loci, IDnum locusCount, char *filename);
void exportContigOrders(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength);
void exportUnusedTranscriptReads(Graph* graph, Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char* directory);
IDnum usedTranscriptReads(Graph * graph, Coordinate minTransLength, Locus * loci, IDnum locusCount);
void exportAMOSTranscripts(Graph * graph, Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char * directory);
void exportTranscriptMappings(Locus * loci, IDnum locusCount, Graph * graph, ReadSet * reads, Coordinate minLength, char * directory);

void removeIndirectConnections();
void cleanTranscriptMemory(Locus * loci, IDnum locusCount);
void cleanLocusMemory(Locus * loci, IDnum locusCount);

void clipTipsHard(Graph * graph, boolean conserveLong);

ReadSet *importEmptyReadSet(char *seqFilename, Coordinate ** lengthsPtr,
			    int wordLength);

void removeNodeFromLocus(Node * node);

void setPairedThreshold(double pairedThreshold);

void exportLocusGraph(FILE * outfile, IDnum index, Locus * loci);

void setDegreeCutoff(int val);

void detachShortReads(ReadSet * reads, int wordLength);
Transcript * allocateTranscript();
#endif
