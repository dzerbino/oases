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
#include <math.h>

#include "globals.h"
#include "recycleBin.h"
#include "utility.h"
#include "graph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "locallyCorrectedGraph.h"
#include "locallyCorrectedGraph2.h"
#include "scaffold.h"
#include "concatenatedGraph.h"
#include "tightString.h"
#include "nodeList.h"
#include "locus.h"
#include "trivialTranscripts.h"
#include "complexTranscripts.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#define BLOCK_SIZE  100000
#define LENGTHCUTOFF 50

struct transcript_st {
	IDnum contigCount;
	Node **contigs;
	Coordinate *distances;
	double confidence;
	Transcript *next;
};

// Global params
static Graph *graph = NULL;

// Global pointers
static RecycleBin *transcriptMemory = NULL;

Transcript *allocateTranscript()
{
	if (transcriptMemory == NULL)
		transcriptMemory =
		    newRecycleBin(sizeof(Transcript), BLOCK_SIZE);

	return allocatePointer(transcriptMemory);
}

Transcript * newTranscript(IDnum contigCount, double confidence) {
	Transcript * transcript = allocateTranscript();
	transcript->contigCount = 0;
	transcript->contigs = callocOrExit(contigCount, Node *);
	transcript->distances = callocOrExit(contigCount, Coordinate);
	transcript->confidence = confidence;
	// DEBUG
	if (confidence > 1)
		abort();
	return transcript;
}

void addContigToTranscript(Transcript * transcript, Node * node, Coordinate distance) {
	transcript->contigs[transcript->contigCount] = node;
	if (transcript->contigCount > 1)
		transcript->distances[transcript->contigCount - 1] = distance;
	transcript->contigCount++;
}

void produceTranscript(Locus * locus, IDnum nodesInList)
{
	IDnum index = 0;
	Node *node;

	Transcript *transcript = newTranscript(nodesInList, ((double) nodesInList) / getContigCount(locus));

	while ((node = popNodeRecord())) {
		transcript->contigs[index] = node;
		if (index > 0) {
			transcript->distances[index - 1] =
			    getConnectionDistance((getConnectionBetweenNodes(transcript->contigs[index - 1], getTwinNode(node))));
			transcript->distances[index - 1] -=
			    getNodeLength(node)/2;
			transcript->distances[index - 1] -=
			    getNodeLength(transcript->contigs[index - 1])/2;
			if (transcript->distances[index - 1] < 0)
				transcript->distances[index - 1] = 0;
		}
		index++;
	}
	transcript->contigCount = index;

	addTranscript(locus, transcript);
}

void setUnreliableConnectionCutoff_oases(int val)
{
	scaffold_setUnreliableConnectionCutoff(val);
}

void cleanTranscriptMemory(Locus * loci, IDnum locusCount)
{
	IDnum index;
	Transcript *transcript;

	for (index = 0; index < locusCount; index++) {
		Locus * locus = getLocus(loci, index);
		for (transcript = getTranscript(locus); transcript;
		     transcript = transcript->next) {
			free(transcript->contigs);
			free(transcript->distances);
		}
	}
	destroyRecycleBin(transcriptMemory);
	transcriptMemory = NULL;
}

static void addTranscriptToLocus(Locus * locus, int configuration, double * scores)
{
	if (configuration == LINEAR)
		addLinearTranscriptToLocus(locus);
	else if (configuration == BUBBLE)
		addBubbleTranscriptsToLocus(locus);
	else if (configuration == FORK)
		addForkedTranscriptsToLocus(locus);
	else if (configuration == UNKNOWN)
		computeHighlyExpressedLocusTranscripts(locus,
						       scores, graph);
}

void computeTranscripts(Graph * argGraph, Locus * loci, IDnum locusCount) {
	IDnum index;
	Locus *locus;
	IDnum distribution[4];
	int configuration;
	double *scores = callocOrExit(2 * nodeCount(argGraph), double);

	graph = argGraph;

	resetNodeStatus(graph);
	prepareGraphForLocalCorrections2(graph);

	for (index = 0; index < locusCount; index++) {
		locus = getLocus(loci, index);
		setLocusStatus(locus, true);
		getDegreeDistribution(locus, distribution);
		configuration = hasPlausibleTranscripts(distribution);
		setLocusStatus(locus, true);
		addTranscriptToLocus(locus, configuration, scores);
		setLocusStatus(locus, false);
	}
	
	deactivateLocalCorrectionSettings2();
	free(scores);
}

void setNextTranscript(Transcript * transcript, Transcript * next) {
	transcript->next = next;
}

Transcript * getNextTranscript(Transcript * transcript) {
	return transcript->next;
}

IDnum getTranscriptContigCount(Transcript * transcript) {
	return transcript->contigCount;
}

Node * getTranscriptContig(Transcript * transcript, IDnum index) {
	return transcript->contigs[index];
}

Coordinate getTranscriptDistance(Transcript * transcript, IDnum index) {
	return transcript->distances[index];
}

double getConfidence(Transcript * transcript) {
	return transcript->confidence;
}
