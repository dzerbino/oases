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

#include "globals.h"
#include "recycleBin.h"
#include "utility.h"
#include "graph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "locallyCorrectedGraph.h"
#include "transcript.h"
#include "scaffold.h"
#include "concatenatedGraph.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#define BLOCK_SIZE  100000
#define LN2 1.4
#define LENGTHCUTOFF 50

typedef enum event_type {
	mutually_exclusive_exons,
	skipped_exon,
	alternative_5prime_splice,
	alternative_3prime_splice,
	intron_retention,
	alternative_polyA,
} EventType;

typedef struct event_st Event;

struct event_st {
	Node *nodes[4];
	EventType type;
	Event *next;
};

struct transcript_st {
	IDnum contigCount;
	Node **contigs;
	Coordinate *distances;
	double confidence;
	Transcript *next;
};

struct locus_st {
	IDnum contigCount;
	IDnum longContigCount;
	Node **contigs;
	Transcript *transcript;
	Event *event;
};

// Global params
static Graph *graph = NULL;
static IDnum maxtrans = 10;
static IDnum UNRELIABLE_CONNECTION_CUTOFF = 5;

// Global pointers
static NodeList *markedNodes;
static RecycleBin *nodeListMemory = NULL;
static RecycleBin *transcriptMemory = NULL;
static RecycleBin *eventMemory = NULL;

// DEBUG
static Node *heavyNode = NULL;

static NodeList *allocateNodeList()
{
	if (nodeListMemory == NULL)
		nodeListMemory =
		    newRecycleBin(sizeof(NodeList), BLOCK_SIZE);

	return allocatePointer(nodeListMemory);
}

static void deallocateNodeList(NodeList * nodeList)
{
	deallocatePointer(nodeListMemory, nodeList);
}

static NodeList *recordNode(Node * node)
{
	NodeList *nodeList = allocateNodeList();
	nodeList->node = node;
	nodeList->next = markedNodes;
	nodeList->previous = NULL;

	if (markedNodes != NULL)
		markedNodes->previous = nodeList;

	markedNodes = nodeList;

	//printf("Recording node %li\n", (long) getNodeID(node));

	return nodeList;
}

static Node *popNodeRecord()
{
	NodeList *nodeList = markedNodes;
	Node *node;

	if (markedNodes == NULL)
		return NULL;

	node = nodeList->node;
	markedNodes = nodeList->next;
	if (markedNodes != NULL)
		markedNodes->previous = NULL;

	//printf("Popping record %li\n", (long) getNodeID(node));

	deallocateNodeList(nodeList);
	return node;
}

static Transcript *allocateTranscript()
{
	if (transcriptMemory == NULL)
		transcriptMemory =
		    newRecycleBin(sizeof(Transcript), BLOCK_SIZE);

	return allocatePointer(transcriptMemory);
}

static void propagateComponent(Node * node)
{
	Connection *connect;

	if (getNodeStatus(node) || !getUniqueness(node))
		return;

	setNodeStatus(node, true);

	for (connect = getConnection(node); connect != NULL;
	     connect = getNextConnection(connect))
		propagateComponent(getConnectionDestination(connect));
	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect))
		propagateComponent(getConnectionDestination(connect));
}

static IDnum countConnectedComponents(Graph * graph)
{
	IDnum index;
	IDnum count = 0;
	Node *node;

	resetNodeStatus(graph);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (!getNodeStatus(node) && getUniqueness(node)) {
			count++;
			propagateComponent(node);
		}
	}

	return count;
}

static void fillUpComponent(Node * node)
{
	Connection *connect;

	if (getNodeStatus(node) || !getUniqueness(node))
		return;
	setSingleNodeStatus(node, true);

	recordNode(node);

	for (connect = getConnection(node); connect != NULL;
	     connect = getNextConnection(connect))
		fillUpComponent(getTwinNode
				(getConnectionDestination(connect)));
	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect))
		fillUpComponent(getConnectionDestination(connect));
}

static IDnum countMarkedNodes()
{
	NodeList *list;
	IDnum counter = 0;

	for (list = markedNodes; list != NULL; list = list->next)
		counter++;

	return counter;
}

static void extendComponentToNode(Node * node)
{
	if (getNodeStatus(node))
		return;

	setSingleNodeStatus(node, true);
	recordNode(node);
}

static void extendComponentFromNode(Node * node)
{
	Connection *connect;
	//printf("Extending from node %li\n", (long) getNodeID(node));

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		extendComponentToNode(getTwinNode
				      (getConnectionDestination(connect)));

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect))
		extendComponentToNode(getConnectionDestination(connect));
}

static void extendComponent(Locus * locus)
{
	IDnum index;

	for (index = 0; index < locus->longContigCount; index++)
		extendComponentFromNode(locus->contigs[index]);
}

static Locus *extractConnectedComponents(IDnum locusCount)
{
	Locus *loci = callocOrExit(locusCount, Locus);
	Locus *locus;
	IDnum index;
	IDnum locusIndex = 0;
	IDnum nodeIndex;
	Node *node;

	resetNodeStatus(graph);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (!getNodeStatus(node) && getUniqueness(node)) {
			nodeIndex = 0;
			locus = &(loci[locusIndex++]);

			// Long contigs
			fillUpComponent(node);
			locus->longContigCount = countMarkedNodes();
			locus->contigs =
			    callocOrExit(locus->longContigCount, Node *);
			while (markedNodes)
				locus->contigs[nodeIndex++] =
				    popNodeRecord();

			// Secondary contigs
			extendComponent(locus);
			locus->contigCount =
			    locus->longContigCount + countMarkedNodes();
			locus->contigs =
			    reallocOrExit(locus->contigs,
					  locus->contigCount, Node *);
			while (markedNodes)
				locus->contigs[nodeIndex++] =
				    popNodeRecord();
			// Mark primary nodes so that their twins are not reused
			for (nodeIndex = 0;
			     nodeIndex < locus->longContigCount;
			     nodeIndex++)
				setNodeStatus(locus->contigs[nodeIndex],
					      true);

			// Unmark secondary nodes so that they are available to other loci
			for (nodeIndex = locus->longContigCount;
			     nodeIndex < locus->contigCount; nodeIndex++)
				setNodeStatus(locus->contigs[nodeIndex],
					      false);
		}
	}

	return loci;
}

static Coordinate isContigInverted(Node * node)
{
	Coordinate nodeLength = getNodeLength(node);
	Coordinate index;
	Nucleotide n1, n2, n3;
	Coordinate totalStops[6];
	// Count of STOP codons per strand/frame
	// 0: + strand, frame 1
	// 1: + strand, frame 2
	// 2: + strand, frame 3
	// 3: - strand, frame 1
	// 4: - strand, frame 2
	// 5: - strand, frame 3
	Coordinate minScorePlus;
	Coordinate minScoreMinus;

	if (getNodeLength(node) < 3)
		return 0;

	// Initialise:
	for (index = 0; index < 6; index++)
		totalStops[index] = 0;
	n2 = getNucleotideInNode(node, 0);
	n3 = getNucleotideInNode(node, 1);

	// Scna sequence and spot Stop codons
	for (index = 2; index < nodeLength; index++) {
		n1 = n2;
		n2 = n3;
		n3 = getNucleotideInNode(node, index);

		if (n1 == THYMINE && n2 == ADENINE
		    && (n3 == ADENINE || n3 == GUANINE))
			totalStops[index % 3]++;

		if (n1 == THYMINE && n2 == GUANINE && n3 == ADENINE)
			totalStops[index % 3]++;

		if (n3 == THYMINE && n2 == ADENINE
		    && (n1 == ADENINE || n1 == GUANINE))
			totalStops[3 + index % 3]++;

		if (n3 == THYMINE && n2 == GUANINE && n1 == ADENINE)
			totalStops[3 + index % 3]++;
	}

	// Find the strand/frame with the minimum number of STOP codons (presumably ORF)
	minScorePlus = totalStops[0];
	for (index = 0; index < 3; index++)
		if (totalStops[index] < minScorePlus)
			minScorePlus = totalStops[index];

	minScoreMinus = totalStops[3];
	for (index = 3; index < 6; index++)
		if (totalStops[index] < minScoreMinus)
			minScoreMinus = totalStops[index];

	if (minScorePlus < minScoreMinus)
		return -nodeLength;
	else if (minScorePlus == minScoreMinus)
		return 0;
	else
		return nodeLength;
}

static boolean isInverted(Locus * locus)
{
	Coordinate score = 0;
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		score += isContigInverted(locus->contigs[index]);

	return score > 0;
}

static void revert(Locus * locus)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		locus->contigs[index] = getTwinNode(locus->contigs[index]);
}

static void orientLoci(Locus * loci, IDnum locusCount)
{
	IDnum index;

	for (index = 0; index < locusCount; index++) {
		if (isInverted(&(loci[index])))
			revert(&(loci[index]));
	}
}

Locus *extractGraphLoci(Graph * argGraph, ReadSet * reads,
			boolean * dubious, Coordinate * lengths,
			IDnum * locusCount)
{
	Locus *loci;

	graph = argGraph;

	buildScaffold(graph, reads, dubious, lengths);

	puts("Extracting loci from connection graph...");

	*locusCount = countConnectedComponents(graph);
	printf("Counted %li mRNA loci\n", (long) *locusCount);

	loci = extractConnectedComponents(*locusCount);
	orientLoci(loci, *locusCount);

	transitiveReduction();

	return loci;
}

static boolean leftHandNeighboursVisited(Node * node)
{
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect))
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect)) ==
		    1)
			return false;

	return true;
}

static IDnum abs_id(IDnum id)
{
	return id > 0 ? id : -id;
}

static Coordinate getTotalCoverage(Node * node)
{
	Category cat;
	Coordinate res = 0;
	PassageMarker *marker;

	for (cat = 0; cat < CATEGORIES; cat++)
		res += getVirtualCoverage(node, cat);

	for (marker = getMarker(node); marker;
	     marker = getNextInNode(marker))
		res += getPassageMarkerLength(marker);

	return res;
}

static void setLocusStatus(Locus * locus, boolean status)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		setSingleNodeStatus(locus->contigs[index], status);
}

static void setNodeConnectionStatus(Node * node, boolean status)
{
	Connection *connect;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			setConnectionStatus(connect, status);
}

static void setLocusConnectionStatus(Locus * locus, boolean status)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		setNodeConnectionStatus(locus->contigs[index], status);
}

static Node *findHeaviestNonUsedNode(Locus * locus)
{
	Node *nextNode = NULL;
	Coordinate maxCov = 0;
	IDnum index;
	Transcript *transcript;

	setLocusStatus(locus, false);

	for (transcript = locus->transcript; transcript;
	     transcript = transcript->next)
		for (index = 0; index < transcript->contigCount; index++)
			setSingleNodeStatus(transcript->contigs[index],
					    true);

	for (index = 0; index < locus->contigCount; index++) {
		if (!getNodeStatus(locus->contigs[index])
		    && getTotalCoverage(locus->contigs[index]) > maxCov) {
			maxCov = getTotalCoverage(locus->contigs[index]);
			nextNode = locus->contigs[index];
		}
	}

	//printf("Next node will be %li\n", (long) getNodeID(nextNode));

	heavyNode = nextNode;
	return nextNode;
}


static void computeDPScore(Node * node, double *scores)
{
	Connection *connect;
	double maxWeight = 0;
	Node *predecessor = NULL;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect))) {
			if (heavyNode && getConnectionDestination(connect) == heavyNode) {
				maxWeight = getConnectionWeight(connect);
				predecessor = getConnectionDestination(connect);
				break;
			} else if (getConnectionWeight(connect) > maxWeight) {
				maxWeight = getConnectionWeight(connect);
				predecessor = getConnectionDestination(connect);
			}
		}
	}

	if (heavyNode && (node == heavyNode || predecessor == heavyNode)) 
		scores[getNodeID(node) + nodeCount(graph)] =
		    100000 * maxWeight + scores[getNodeID(predecessor) + nodeCount(graph)];
	else
		scores[getNodeID(node) + nodeCount(graph)] =
		    maxWeight + scores[getNodeID(predecessor) + nodeCount(graph)];
}

static void propagateDP(Node * node)
{
	Connection *connect;
	Node *nextNode;

	for (connect = getConnection(node); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (!getConnectionStatus(connect))
			continue;

		nextNode = getTwinNode(getConnectionDestination(connect));
		if (getNodeStatus(nextNode) != 1)
			continue;

		if (leftHandNeighboursVisited(nextNode))
			recordNode(nextNode);
	}
}

static Node *chooseBestUnvisitedPredecessor(Node * node)
{
	Node *nextNode = NULL;
	double maxWeight = 0;
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect)) ==
		    1 && getConnectionWeight(connect) > maxWeight) {
			nextNode = getConnectionDestination(connect);
			maxWeight = getConnectionWeight(connect);
		}
	}

	return nextNode;
}

static Connection *getReverseActiveConnection(Node * node)
{
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus(getConnectionDestination(connect)))
			return connect;

	return false;
}

static Connection *getReverseMarkedConnection(Node * node)
{
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect))
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect)))
			return connect;

	return false;
}

static void destroyNodeBackConnections(Node * node)
{
	Connection *connect, *next;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = next) {
		next = getNextConnection(connect);
		if (getNodeStatus(getConnectionDestination(connect)) == 3)
			setConnectionStatus(connect, false);
	}
}

static void createAlternativeDPStartFromNode(Node * node)
{
	Node *current = node;
	Node *next;

	setSingleNodeStatus(current, 3);

	while ((next = chooseBestUnvisitedPredecessor(current))) {
		current = next;
		setSingleNodeStatus(current, 3);
	}

	if (getConnectionBetweenNodes(node, getTwinNode(current)))
		destroyNodeBackConnections(node);
	else
		destroyNodeBackConnections(current);
}

static boolean lookForPossibleAlternativeDPStartFromNode(Node * node)
{
	Connection *connect;

	if (getNodeStatus(node) != 1)
		return false;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect)) ==
		    2) {
			createAlternativeDPStartFromNode(node);
			return true;
		}
	}

	return false;
}

static boolean createAlternativeDPStarts(Locus * locus)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		if (lookForPossibleAlternativeDPStartFromNode
		    (locus->contigs[index]))
			return true;

	return false;
}

static Node *chooseBestPredecessor(Node * node)
{
	Node *nextNode = NULL;
	double maxWeight = 0;
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (!getConnectionStatus(connect))
			continue;

		if (heavyNode 
		    && getConnectionDestination(connect) == heavyNode
		    && getNodeStatus(heavyNode) == 2)
			return heavyNode;

		if (getNodeStatus(getConnectionDestination(connect)) == 2
		    && getConnectionWeight(connect) > maxWeight) {
			nextNode = getConnectionDestination(connect);
			maxWeight = getConnectionWeight(connect);
		}
	}

	return nextNode;
}

static IDnum extractMajorityPath(Node * maxNode)
{
	Node *node = maxNode;
	IDnum nodesInList = 1;

	recordNode(node);
	setSingleNodeStatus(node, 1);

	while (getReverseMarkedConnection(node)) {
		//printf("Back tracking through node %li\n", (signed long) getNodeID(node));
		node = chooseBestPredecessor(node);
		if (node == NULL)
			break;
		recordNode(node);
		setSingleNodeStatus(node, 1);
		nodesInList++;
	}
	//printf("Arrived in %li\n", (long) getNodeID(node));

	return nodesInList;
}

static IDnum extractLinearPath(Node * maxNode)
{
	Node *node = maxNode;
	IDnum nodesInList = 1;
	Connection *connect;

	recordNode(node);
	setNodeStatus(node, false);

	while ((connect = getReverseActiveConnection(node))) {
		//printf("Back tracking through node %li\n", (signed long) getNodeID(node));
		node = getConnectionDestination(connect);
		recordNode(node);
		setNodeStatus(node, false);
		nodesInList++;
	}
	//printf("Arrived in %li\n", (long) getNodeID(node));

	return nodesInList;
}

static void produceTranscript(Locus * locus, IDnum nodesInList)
{
	IDnum index = 0;
	Node *node;

	Transcript *transcript = allocateTranscript();
	transcript->contigCount = nodesInList;
	transcript->contigs = callocOrExit(nodesInList, Node *);
	transcript->distances = callocOrExit(nodesInList, Coordinate);
	transcript->confidence = 1;

	while ((node = popNodeRecord())) {
		transcript->contigs[index] = node;
		if (index > 0) {
			transcript->distances[index - 1] =
			    getConnectionDistance((getConnectionBetweenNodes(transcript->contigs[index - 1], getTwinNode(node))));
			transcript->distances[index - 1] -=
			    getNodeLength(node);
			transcript->distances[index - 1] -=
			    getNodeLength(transcript->contigs[index - 1]);
			if (transcript->distances[index - 1] < 0)
				transcript->distances[index - 1] = 0;
		}
		index++;
	}
	transcript->next = locus->transcript;
	locus->transcript = transcript;
}

void computeHighestExpressedLocusTranscript(Locus * locus, double *scores)
{
	IDnum index;
	Node *node;
	double overallMaxScore = -1;
	Node *maxNode = NULL;
	IDnum nodesToVisit = locus->contigCount;
	IDnum nodesInList;

	// Clear node statuses:
	// Record all nodes from which to start the DP
	setLocusStatus(locus, true);

	for (index = 0; index < locus->contigCount; index++) {
		node = locus->contigs[index];
		if (getReverseMarkedConnection(node) == NULL)
			recordNode(node);
	}

	// OK, let's name a volunteer...
	if (markedNodes == NULL) {
		createAlternativeDPStartFromNode(locus->contigs[0]);
		computeHighestExpressedLocusTranscript(locus, scores);
		return;
	}
	// Propagate DP
	while ((node = popNodeRecord())) {
		//printf("Visiting node %li\n", (long) getNodeID(node));
		setSingleNodeStatus(node, 2);
		nodesToVisit--;

		computeDPScore(node, scores);
		if (scores[getNodeID(node) + nodeCount(graph)] > overallMaxScore
		    && (!heavyNode || node != getTwinNode(heavyNode))) {
			overallMaxScore = scores[getNodeID(node) + nodeCount(graph)];
			maxNode = node;
		}
		propagateDP(node);

		// If loop prevents nodes from being visited:
		if (nodesToVisit > 0 && markedNodes == NULL) {
			while (popNodeRecord());
			if (createAlternativeDPStarts(locus))
				computeHighestExpressedLocusTranscript
				    (locus, scores);
			return;
		}
	}

	nodesInList = extractMajorityPath(maxNode);
	produceTranscript(locus, nodesInList);
	locus->transcript->confidence =
	    locus->transcript->contigCount / (double) locus->contigCount;
	findHeaviestNonUsedNode(locus);
}

static IDnum explainedContigs(Locus * locus)
{
	IDnum index;
	Node *node;
	Transcript *transcript;
	IDnum counter = 0;

	setLocusStatus(locus, false);

	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next) {
		for (index = 0; index < transcript->contigCount; index++) {
			node = transcript->contigs[index];
			if (!getNodeStatus(node)) {
				counter++;
				setSingleNodeStatus(node, true);
			}
		}
	}

	//printf("%li nodes explained out of %li\n", (long) counter, (long) locus->contigCount);

	return counter;
}

static void computeHighlyExpressedLocusTranscripts(Locus * locus,
						   double *scores)
{
	IDnum counter = 0;

	heavyNode = NULL;

	setLocusStatus(locus, true);
	setLocusConnectionStatus(locus, true);

	while (counter++ < maxtrans
	       && explainedContigs(locus) < locus->contigCount)
		computeHighestExpressedLocusTranscript(locus, scores);


	setLocusStatus(locus, true);
	setLocusConnectionStatus(locus, false);
	setLocusStatus(locus, false);
}

void computeHighlyExpressedTranscripts(Locus * loci, IDnum locusCount)
{
	IDnum index;
	double *scores = callocOrExit(2 * nodeCount(graph), double);

	puts("Computing highly expressed loci...");

	resetNodeStatus(graph);

	for (index = 0; index < locusCount; index++)
		computeHighlyExpressedLocusTranscripts(&(loci[index]),
						       scores);

	free(scores);
}

void setUnreliableConnectionCutoff_oases(int val)
{
	UNRELIABLE_CONNECTION_CUTOFF = (IDnum) val;
}

void cleanTranscriptMemory(Locus * loci, IDnum locusCount)
{
	IDnum index;
	Transcript *transcript;

	for (index = 0; index < locusCount; index++) {
		for (transcript = loci[index].transcript; transcript;
		     transcript = transcript->next) {
			free(transcript->contigs);
			free(transcript->distances);
		}
		loci[index].transcript = NULL;
	}

}

void cleanLocusMemory(Locus * loci, IDnum locusCount)
{
	IDnum index;

	for (index = 0; index < locusCount; index++)
		free(loci[index].contigs);
	free(loci);
	destroyRecycleBin(transcriptMemory);
	destroyRecycleBin(nodeListMemory);
	destroyRecycleBin(eventMemory);
	nodeListMemory = NULL;
	transcriptMemory = NULL;
	eventMemory = NULL;
	cleanScaffoldMemory();
}

static void exportContigSequence(Node * node, FILE * outfile, int *column,
				 boolean hideKmer)
{
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	Coordinate start = 0;
	char str[100];

	if (hideKmer && getNodeLength(node) >= getWordLength(graph) - 1)
		start = getWordLength(graph) - 1;

	while (start < strlen(string)) {
		if (getNodeLength(node) < getWordLength(graph) - 1)
			strncpy(str, string + start, 60 - *column);
		else
			strncpy(str, string + start, 60 - *column);
		str[60 - *column] = '\0';
		fprintf(outfile, "%s", str);
		*column += strlen(str);
		if (*column >= 60) {
			fprintf(outfile, "\n");
			*column = 0;
		}
		start += strlen(str);
	}

	free(string);
}

static void exportGapSequence(Coordinate length, FILE * outfile,
			      int *column)
{
	IDnum index;

	for (index = 0; index < length; index++) {
		fprintf(outfile, "N");

		if ((*column)++ == 60) {
			fprintf(outfile, "\n");
			*column = 0;
		}
	}
}

static void exportTranscript(Transcript * transcript, IDnum locusID,
			     IDnum transID, FILE * outfile)
{
	IDnum index;
	int column = 0;
	boolean hideKmer = false;

	// Header
	fprintf(outfile, ">Locus_%li_Transcript_%li_Confidence_%.3f\n",
		(long) locusID, (long) transID, transcript->confidence);

	// Sequence
	for (index = 0; index < transcript->contigCount; index++) {
		exportContigSequence(transcript->contigs[index], outfile,
				     &column, hideKmer);
		if (index < transcript->contigCount - 1) {
			exportGapSequence(transcript->distances[index],
					  outfile, &column);
			if (transcript->distances[index] <= 0)
				hideKmer = true;
			else
				hideKmer = false;
		}

	}

	if (column > 0)
		fprintf(outfile, "\n");
}

static void exportLocusTranscripts(Locus * locus, IDnum locusID,
				   FILE * outfile)
{
	IDnum index = 0;
	Transcript *transcript;

	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next)
		exportTranscript(transcript, locusID, index++, outfile);
}

void exportTranscripts(Locus * loci, IDnum locusCount, char *filename)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;

	printf("Exporting transcripts to %s\n", filename);

	for (index = 0; index < locusCount; index++)
		exportLocusTranscripts(&(loci[index]), index, outfile);
	fclose(outfile);
}

static IDnum getDegree(Node * node)
{
	Connection *connect;
	IDnum counter = 0;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			counter++;

	return counter;
}

static IDnum getReverseDegree(Node * node)
{
	Connection *connect;
	IDnum counter = 0;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus(getConnectionDestination(connect)))
			counter++;

	return counter;
}

static void updateDistribution(IDnum * distribution, Node * node)
{
	IDnum degree = getDegree(node);

	if (degree < 4)
		distribution[degree]++;
	else
		distribution[3]++;

	degree = getReverseDegree(getTwinNode(node));

	if (degree < 4)
		distribution[degree]++;
	else
		distribution[3]++;
}

static void getDegreeDistribution(Locus * locus, IDnum * distribution)
{
	IDnum index;

	for (index = 0; index < 4; index++)
		distribution[index] = 0;

	for (index = 0; index < locus->contigCount; index++)
		updateDistribution(distribution, locus->contigs[index]);
}

#define UNKNOWN 0
#define LINEAR 1
#define FORK 2
#define BUBBLE 3

static boolean hasPlausibleTranscripts(IDnum * d)
{
	// Singleton or linear chain
	if (d[0] == 2 && d[2] == 0 && d[3] == 0)
		return LINEAR;

	// Alternative start / end
	if (d[0] == 3 && d[2] == 1 && d[3] == 0)
		return FORK;

	// Skipped exon/ 5' or 3' alternative site/ intron retention 
	if (d[0] == 2 && d[2] == 2 && d[3] == 0)
		return BUBBLE;

	return UNKNOWN;
}

static boolean hasNoActiveConnections(Node * node)
{
	Connection *connect;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			return false;

	return true;
}

static Node *findEndNodeInLocus(Locus * locus)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		if (hasNoActiveConnections(locus->contigs[index]))
			return locus->contigs[index];

	return NULL;
}

static void addLinearTranscriptToLocus(Locus * locus)
{
	Node *node = findEndNodeInLocus(locus);
	IDnum nodesInList = extractLinearPath(node);
	produceTranscript(locus, nodesInList);
}

static Connection *getActiveConnection(Node * node)
{
	Connection *connect;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			return connect;

	return false;
}

static Connection *getSecondActiveConnection(Node * node)
{
	Connection *connect;
	boolean firstFound = false;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect)) {
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect)))) {
			if (firstFound)
				return connect;
			else
				firstFound = true;
		}
	}

	return false;
}

static Connection *getSecondReverseActiveConnection(Node * node)
{
	Connection *connect;
	boolean firstFound = false;

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect)) {
		if (getNodeStatus(getConnectionDestination(connect))) {
			if (firstFound)
				return connect;
			else
				firstFound = true;
		}
	}

	return false;
}

static IDnum countActiveConnections(Node * node)
{
	Connection *connect;
	IDnum counter = 0;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			counter++;

	return counter;
}

static void addBubbleTranscriptsToLocus(Locus * locus)
{
	Node *node = findEndNodeInLocus(locus);
	IDnum count = 0;

	// Transcript 1:
	while (node) {
		recordNode(node);
		setSingleNodeStatus(node, false);
		if (getReverseActiveConnection(node))
			node =
			    getConnectionDestination
			    (getReverseActiveConnection(node));
		else
			node = NULL;
		count++;
	}

	produceTranscript(locus, count);
	setLocusStatus(locus, true);

	// Transcript 2:
	node = findEndNodeInLocus(locus);
	count = 0;
	while (node) {
		recordNode(node);
		setSingleNodeStatus(node, false);
		if (getSecondReverseActiveConnection(node))
			node =
			    getConnectionDestination
			    (getSecondReverseActiveConnection(node));
		else if (getReverseActiveConnection(node))
			node =
			    getConnectionDestination
			    (getReverseActiveConnection(node));
		else
			node = NULL;
		count++;
	}

	produceTranscript(locus, count);
}

static Node *findSecondEndNodeInLocus(Locus * locus)
{
	IDnum index;
	boolean firstFound = false;

	for (index = 0; index < locus->contigCount; index++) {
		if (hasNoActiveConnections(locus->contigs[index])) {
			if (firstFound)
				return locus->contigs[index];
			else
				firstFound = true;
		}
	}

	return NULL;
}

static void add3primeForkedTranscriptsToLocus(Locus * locus)
{
	Node *node = findEndNodeInLocus(locus);
	IDnum nodesInList = extractLinearPath(node);
	produceTranscript(locus, nodesInList);

	setLocusStatus(locus, true);

	node = findSecondEndNodeInLocus(locus);
	nodesInList = extractLinearPath(node);
	produceTranscript(locus, nodesInList);
}

static boolean onlyOneEndPoint(Locus * locus)
{
	return findSecondEndNodeInLocus(locus) == NULL;
}

static void addForkedTranscriptsToLocus(Locus * locus)
{
	if (onlyOneEndPoint(locus))
		addBubbleTranscriptsToLocus(locus);
	else
		add3primeForkedTranscriptsToLocus(locus);
}

static void addPlausibleTranscriptToLocus(Locus * locus, int configuration)
{
	if (configuration == LINEAR)
		addLinearTranscriptToLocus(locus);
	else if (configuration == BUBBLE)
		addBubbleTranscriptsToLocus(locus);
	else if (configuration == FORK)
		addForkedTranscriptsToLocus(locus);
	else if (configuration == UNKNOWN);
}

void computePlausibleTranscripts(Locus * loci, IDnum locusCount)
{
	IDnum index;
	Locus *locus;
	IDnum distribution[4];
	int configuration;

	resetNodeStatus(graph);

	for (index = 0; index < locusCount; index++) {
		locus = &(loci[index]);
		setLocusStatus(locus, true);
		getDegreeDistribution(locus, distribution);
		configuration = hasPlausibleTranscripts(distribution);
		addPlausibleTranscriptToLocus(locus, configuration);
		setLocusStatus(locus, false);
	}
}

static void removeIndirectConnectionsAtIndex(IDnum index)
{
	Connection *connect = getConnection(getNodeInGraph(graph, index));
	Connection *next;

	while (connect) {
		next = getNextConnection(connect);
		if (getConnectionDirectCount(connect) == 0)
			destroyConnection(connect, index);
		connect = next;
	}
}

void removeIndirectConnections()
{
	IDnum index;

	for (index = -nodeCount(graph); index <= nodeCount(graph); index++)
		removeIndirectConnectionsAtIndex(index);
}

static Event *allocateEvent()
{
	if (eventMemory == NULL)
		eventMemory = newRecycleBin(sizeof(Event), BLOCK_SIZE);

	return allocatePointer(eventMemory);
}

#define SPLICE_FUZZINESS 3

static boolean donorSiteAtJunction(Node * nodeA, Node * nodeB)
{
	Nucleotide n1, n2;
	int i;

	n2 = getNucleotideInNode(nodeA,
				 getNodeLength(nodeA) - SPLICE_FUZZINESS);

	for (i = SPLICE_FUZZINESS - 1; i > 0; i--) {
		n1 = n2;
		n2 = getNucleotideInNode(nodeA, getNodeLength(nodeA) - i);
		if (n1 == GUANINE && n2 == THYMINE)
			return true;
	}

	for (i = 0; i < SPLICE_FUZZINESS + 2; i++) {
		n1 = n2;
		n2 = getNucleotideInNode(nodeB, i);
		if (n1 == GUANINE && n2 == THYMINE)
			return true;
	}

	return false;
}

static boolean acceptorSiteAtJunction(Node * nodeA, Node * nodeB)
{
	Node *twinNodeA = getTwinNode(nodeA);
	Node *twinNodeB = getTwinNode(nodeB);
	Nucleotide n1, n2;
	int i;

	n2 = getNucleotideInNode(twinNodeB,
				 getNodeLength(twinNodeB) -
				 SPLICE_FUZZINESS);

	for (i = SPLICE_FUZZINESS - 1; i > 0; i--) {
		n1 = n2;
		n2 = getNucleotideInNode(twinNodeB,
					 getNodeLength(twinNodeB) - i);
		if (n1 == CYTOSINE && n2 == ADENINE)
			return true;
	}

	for (i = 0; i < SPLICE_FUZZINESS + 2; i++) {
		n1 = n2;
		n2 = getNucleotideInNode(twinNodeA, i);
		if (n1 == CYTOSINE && n2 == ADENINE)
			return true;
	}

	return false;
}

static boolean finishesWithPAS(Node * node)
{
	char *nodeSeq = expandNodeFragment(node, 0, getNodeLength(node),
					   getWordLength(graph));
	boolean res = false;

	char *ptr = strstr(nodeSeq, "AATAAA");
	if (ptr)
		res = true;
	ptr = strstr(nodeSeq, "ATTAAA");
	if (ptr)
		res = true;

	free(nodeSeq);
	return res;
}

static void extractNodeASEvents(Node * node, Locus * locus)
{
	Node *nodeA, *nodeB, *nodeC;
	Event *event;

	// If linear or more than 2 outgoing arcs: ignore
	if (countActiveConnections(node) != 2)
		return;

	// Follow the two active arcs
	nodeA =
	    getTwinNode(getConnectionDestination
			(getActiveConnection(node)));
	nodeB =
	    getTwinNode(getConnectionDestination
			(getSecondActiveConnection(node)));

	// A should be the longer of the two
	if (getNodeLength(nodeA) < getNodeLength(nodeB)) {
		nodeC = nodeA;
		nodeA = nodeB;
		nodeB = nodeC;
		nodeC = NULL;
	}
	// If both very short, ignore:
	if (getNodeLength(nodeA) < 2 * getWordLength(graph) - 1)
		return;

	if (getNodeLength(nodeB) < 2 * getWordLength(graph) - 1) {
		if (countActiveConnections(nodeA) != 1
		    || countActiveConnections(nodeB) != 1
		    || getConnectionDestination(getActiveConnection(nodeA))
		    !=
		    getConnectionDestination(getActiveConnection(nodeB)))
			return;

		nodeC =
		    getTwinNode(getConnectionDestination
				(getActiveConnection(nodeA)));

		// Intron retention
		if (donorSiteAtJunction(node, nodeA)
		    && acceptorSiteAtJunction(nodeA, nodeC)) {
			event = allocateEvent();
			event->type = intron_retention;
			event->nodes[0] = node;
			event->nodes[1] = nodeA;
			event->nodes[2] = nodeB;
			event->nodes[3] = nodeC;
			event->next = locus->event;
			locus->event = event;
		}
		// Alternative 5' splice site
		else if (donorSiteAtJunction(node, nodeA)) {
			event = allocateEvent();
			event->type = alternative_5prime_splice;
			event->nodes[0] = node;
			event->nodes[1] = nodeA;
			event->nodes[2] = nodeB;
			event->nodes[3] = nodeC;
			event->next = locus->event;
			locus->event = event;
		}
		// Alternative 3' splice site
		else if (acceptorSiteAtJunction(nodeA, nodeC)) {
			event = allocateEvent();
			event->type = alternative_3prime_splice;
			event->nodes[0] = node;
			event->nodes[1] = nodeA;
			event->nodes[2] = nodeB;
			event->nodes[3] = nodeC;
			event->next = locus->event;
			locus->event = event;
		}
		// Skipped exon
		else {
			event = allocateEvent();
			event->type = skipped_exon;
			event->nodes[0] = node;
			event->nodes[1] = nodeA;
			event->nodes[2] = nodeB;
			event->nodes[3] = nodeC;
			event->next = locus->event;
			locus->event = event;
		}
	} else {
		// Alt. poly A:
		if (finishesWithPAS(node) && finishesWithPAS(nodeA)) {
			event = allocateEvent();
			event->type = alternative_polyA;
			event->nodes[0] = node;
			event->nodes[1] = nodeA;
			event->nodes[2] = nodeB;
			event->nodes[3] = NULL;
			event->next = locus->event;
			locus->event = event;
		}
		// Mutually exclusive exons
		if (countActiveConnections(nodeA) == 1
		    && countActiveConnections(nodeB) == 1
		    && getConnectionDestination(getActiveConnection(nodeA))
		    ==
		    getConnectionDestination(getActiveConnection(nodeB))) {
			event = allocateEvent();
			event->type = mutually_exclusive_exons;
			event->nodes[0] = node;
			event->nodes[1] = nodeA;
			event->nodes[2] = nodeB;
			event->nodes[3] =
			    getTwinNode(getConnectionDestination
					(getActiveConnection(nodeA)));
			event->next = locus->event;
			locus->event = event;
		}
	}
}

static void extractLocusASEvents(Locus * locus)
{
	IDnum index;

	setLocusStatus(locus, true);

	for (index = 0; index < locus->longContigCount; index++)
		extractNodeASEvents(locus->contigs[index], locus);

	setLocusStatus(locus, false);
}

void computeASEvents(Locus * loci, IDnum locusCount)
{
	IDnum index;

	puts("Extracting identifiable AS events");

	resetNodeStatus(graph);
	for (index = 0; index < locusCount; index++)
		extractLocusASEvents(&(loci[index]));
}

static void exportNodeSequence(Node * node, FILE * outfile)
{
	int wordShift = getWordLength(graph) - 1;
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	Coordinate length = getNodeLength(node) + wordShift;
	Coordinate start;
	char str[100];

	for (start = 0; start <= length; start += 60) {
		strncpy(str, &(string[start]), 60);
		str[60] = '\0';
		fprintf(outfile, "%s\n", str);
	}
	free(string);
}

static void exportLocusNode(IDnum index, Node * node, FILE * outfile)
{
	fprintf(outfile, ">Locus_%li_Node_%li\n", (long) index,
		(long) abs_id(getNodeID(node)));
	exportNodeSequence(node, outfile);
}

static void exportLocusNodes(IDnum locusIndex, Locus * locus,
			     FILE * outfile)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		exportLocusNode(locusIndex, locus->contigs[index],
				outfile);
}

static void getNodeStringID(char *str, Node * node)
{
	IDnum id = getNodeID(node);

	sprintf(str, "%li", (long) abs_id(id));
}

static void exportEvent(IDnum index, Event * event, FILE * outfile)
{
	char id[4][100];
	int i;

	for (i = 0; i < 4; i++)
		getNodeStringID(id[i], event->nodes[i]);

	fprintf(outfile, "Locus %li: ", (long) index);

	if (event->type == mutually_exclusive_exons) {
		fprintf(outfile,
			"[MEE] Mutually exclusive exons : %sS-%sE^ (%sS-%sE^,%sS-%sE^) %sS-%sE^\n",
			id[0],
			id[0], id[1], id[1], id[2], id[2], id[3], id[3]);
	} else if (event->type == skipped_exon) {
		fprintf(outfile,
			"[SE] Skipped exon : %sS-%sE^ (%sS-%sE^,0) %sS-%sE^\n",
			id[0], id[0], id[1], id[1], id[3], id[3]);
	} else if (event->type == alternative_5prime_splice) {
		fprintf(outfile,
			"[a5p] Alternative 5' splice: %sS- (%sE^, %sE^) %sS-%sE^\n",
			id[0], id[0], id[1], id[3], id[3]);
	} else if (event->type == alternative_3prime_splice) {
		fprintf(outfile,
			"[a3p] Alternative 3' splice: %sS-%sE^ (%sS-,%sS-) %sE^\n",
			id[0], id[0], id[1], id[3], id[3]);
	} else if (event->type == intron_retention) {
		fprintf(outfile,
			"[IR] Intron retention: %sS- (%sS^%sS-,0) %sE^\n",
			id[0], id[0], id[3], id[3]);
	} else if (event->type == alternative_polyA) {
		fprintf(outfile,
			"[aPS] Alternative polyadenylation site: %sS (-%sE], -%sE])\n",
			id[0], id[0], id[1]);
	}
}

static void exportLocusASEvents(IDnum index, Locus * locus, FILE * outfile)
{
	Event *event;

	exportLocusNodes(index, locus, outfile);
	for (event = locus->event; event; event = event->next)
		exportEvent(index, event, outfile);
}

void exportASEvents(Locus * loci, IDnum locusCount, char *filename)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;

	if (outfile)
		printf("Exporting AS events to %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	for (index = 0; index < locusCount; index++)
		exportLocusASEvents(index, &(loci[index]), outfile);

	fclose(outfile);
}

static void exportTranscriptContigs(Transcript * transcript, IDnum locusID,
				    IDnum transID, FILE * outfile)
{
	IDnum index;
	Coordinate totalLength = getWordLength(graph) - 1;

	// Header
	fprintf(outfile, ">Locus_%li_Transcript_%li_Confidence_%.3f\n",
		(long) locusID, (long) transID, transcript->confidence);

	// Sequence
	for (index = 0; index < transcript->contigCount; index++) {
		totalLength += getNodeLength(transcript->contigs[index]);
		fprintf(outfile, "%li:%lli",
			(long) getNodeID(transcript->contigs[index]),
			(long long) totalLength);
		if (index < transcript->contigCount - 1) {
			fprintf(outfile, "-(%li)->",
				(long) transcript->distances[index]);
			totalLength += transcript->distances[index];
		}
	}

	fprintf(outfile, "\n");
}

static void exportLocusContigs(IDnum locusID, Locus * locus,
			       FILE * outfile)
{
	IDnum index = 0;
	Transcript *transcript;

	exportLocusNodes(locusID, locus, outfile);
	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next)
		exportTranscriptContigs(transcript, locusID, index++,
					outfile);
}

void exportContigOrders(Locus * loci, IDnum locusCount, char *filename)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;

	if (outfile)
		printf("Exporting transcript contigs to %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	for (index = 0; index < locusCount; index++)
		exportLocusContigs(index, &(loci[index]), outfile);

	fclose(outfile);

}

void exportLocusGraph(FILE * file, IDnum index, Locus * loci)
{
	Locus *locus = &loci[index];
	IDnum i;
	long long nodeID;
	Connection *connect;
	Node *node;

	resetNodeStatus(graph);
	setLocusStatus(locus, true);

	fprintf(file, "digraph graph_%lli {\n", (long long) index);

	fprintf(file, "rankdir=LR\n");
	fprintf(file, "style=invis\n");
	fprintf(file, "node [shape=point]\n");
	fprintf(file, "edge [style=bold, color=red]\n");
	for (i = 0; i < locus->longContigCount; i++) {
		nodeID = (long long) getNodeID(locus->contigs[i]);
		if (nodeID > 0)
			fprintf(file, "subgraph cluster%li {\n",
				(long) nodeID);
		else
			fprintf(file, "subgraph cluster0%li {\n",
				(long) -nodeID);
		fprintf(file, "%lli -> %lli [label=%lli]\n", nodeID,
			-nodeID, nodeID);
		fprintf(file, "}\n");
	}

	fprintf(file, "edge [color=black]\n");
	for (i = locus->longContigCount; i < locus->contigCount; i++) {
		nodeID = (long long) getNodeID(locus->contigs[i]);
		if (nodeID > 0)
			fprintf(file, "subgraph cluster%li {\n",
				(long) nodeID);
		else
			fprintf(file, "subgraph cluster0%li {\n",
				(long) -nodeID);
		fprintf(file, "%lli -> %lli [label=%lli]\n", nodeID,
			-nodeID, nodeID);
		fprintf(file, "}\n");
	}

	fprintf(file, "edge [style=normal]\n");
	for (i = 0; i < locus->contigCount; i++) {
		node = locus->contigs[i];
		nodeID = getNodeID(node);

		for (connect = getConnection(node); connect != NULL;
		     connect = getNextConnection(connect))
			if (getNodeStatus
			    (getTwinNode
			     (getConnectionDestination(connect)))) {
				fprintf(file, "%lli -> %lli [label=%lli",
					-nodeID,
					(long long)
					-getNodeID(getConnectionDestination
						   (connect)),
					(long long)
					getConnectionDistance(connect));
				if (!getConnectionStatus(connect))
					fprintf(file, "XXX");
				fprintf(file, "]\n");
			}
	}

	fprintf(file, "}\n");

	setLocusStatus(locus, false);
}

ReadSet *importEmptyReadSet(char *filename, Coordinate ** lengthsPtr,
			    int wordLength)
{
	FILE *file = fopen(filename, "r");
	const int maxline = 5000;
	char line[5000];
	IDnum sequenceCount, sequenceIndex;
	ReadSet *reads;
	short int temp_short;
	int lengthOffset = wordLength - 1;
	Coordinate bpCount = 0;

	if (file != NULL)
		printf("Reading read set file %s;\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	reads = newReadSet();

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);
	printf("%d sequences found\n", sequenceCount);

	reads->readCount = sequenceCount;

	if (reads->readCount == 0) {
		reads->sequences = NULL;
		reads->categories = NULL;
		return reads;
	}

	*lengthsPtr = callocOrExit(sequenceCount, Coordinate);
	reads->categories = callocOrExit(sequenceCount, Category);
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {

			// Reading category info
			sscanf(line, "%*[^\t]\t%*[^\t]\t%hd", &temp_short);
			reads->categories[sequenceIndex + 1] =
			    (Category) temp_short;
			if (sequenceIndex != -1)
				(*lengthsPtr)[sequenceIndex] =
				    bpCount - lengthOffset;
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (Coordinate) strlen(line) - 1;
		}
	}
	(*lengthsPtr)[sequenceIndex] = bpCount - lengthOffset;
	fclose(file);

	puts("Done");
	return reads;

}

boolean *removeLowCoverageNodesAndDenounceDubiousReads(Graph * graph,
						       double minCov)
{
	IDnum index;
	Node *node;
	boolean denounceReads = readStartsAreActivated(graph);
	boolean *res = NULL;
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	IDnum maxIndex;
	IDnum readID;
	IDnum index2;

	printf("Removing contigs with coverage < %f...\n", minCov);

	if (denounceReads)
		res = callocOrExit(sequenceCount(graph), boolean);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);

		if (getNodeLength(node) == 0)
			continue;

		if (getTotalCoverage(node) / getNodeLength(node) < minCov) {
			if (denounceReads) {
				nodeArray = getNodeReads(node, graph);
				maxIndex = getNodeReadCount(node, graph);
				for (index2 = 0; index2 < maxIndex;
				     index2++) {
					shortMarker =
					    getShortReadMarkerAtIndex
					    (nodeArray, index2);
					readID =
					    getShortReadMarkerID
					    (shortMarker);
					//printf("Dubious %d\n", readID);
					if (readID > 0)
						res[readID - 1] = true;
					else
						res[-readID - 1] = true;
				}

				nodeArray =
				    getNodeReads(getTwinNode(node), graph);
				maxIndex =
				    getNodeReadCount(getTwinNode(node),
						     graph);
				for (index2 = 0; index2 < maxIndex;
				     index2++) {
					shortMarker =
					    getShortReadMarkerAtIndex
					    (nodeArray, index2);
					readID =
					    getShortReadMarkerID
					    (shortMarker);
					//printf("Dubious %d\n", readID);
					if (readID > 0)
						res[readID - 1] = true;
					else
						res[-readID - 1] = true;
				}
			}

			while ((marker = getMarker(node))) {
				if (!isInitial(marker)
				    && !isTerminal(marker))
					disconnectNextPassageMarker
					    (getPreviousInSequence(marker),
					     graph);
				destroyPassageMarker(marker);
			}
			destroyNode(node, graph);
		}
	}

	concatenateGraph(graph);
	return res;
}

static Coordinate getTipLength(Node * node)
{
	Node *current = getTwinNode(node);
	Coordinate length = 0;

	if (simpleArcCount(current) > 1)
		return getNodeLength(node);

	while (current != NULL && simpleArcCount(getTwinNode(current)) < 2
	       && simpleArcCount(current) < 2) {
		length += getNodeLength(current);
		current = getDestination(getArc(current));
	}

	return length;
}

void clipTipsHard(Graph * graph)
{
	IDnum index;
	Node *current, *twin;
	boolean modified = true;
	int Wordlength = getWordLength(graph);
	PassageMarker *marker;

	puts("Clipping short tips off graph, drastic");

	while (modified) {
		modified = false;
		for (index = 1; index <= nodeCount(graph); index++) {
			current = getNodeInGraph(graph, index);

			if (current == NULL)
				continue;

			twin = getTwinNode(current);

			if (getArc(current) == NULL
			    && getTipLength(current) < 2 * Wordlength) {
				while ((marker = getMarker(current))) {
					if (!isInitial(marker)
					    && !isTerminal(marker))
						disconnectNextPassageMarker
						    (getPreviousInSequence
						     (marker), graph);
					destroyPassageMarker(marker);
				}
				destroyNode(current, graph);
				modified = true;
			} else if (getArc(twin) == NULL
				   && getTipLength(twin) <
				   2 * Wordlength) {
				while ((marker = getMarker(current))) {
					if (!isInitial(marker)
					    && !isTerminal(marker))
						disconnectNextPassageMarker
						    (getPreviousInSequence
						     (marker), graph);
					destroyPassageMarker(marker);
				}
				destroyNode(twin, graph);
				modified = true;
			}
		}
	}

	concatenateGraph(graph);
	printf("%d nodes left\n", nodeCount(graph));
}
