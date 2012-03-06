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

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#define BLOCK_SIZE  100000
#define LENGTHCUTOFF 50

typedef enum event_type {
	mutually_exclusive_exons,
	skipped_exon,
	alternative_5prime_splice,
	alternative_3prime_splice,
	intron_retention,
	alternative_polyA,
} EventType;

struct event_st {
	Node *nodes[4];
	EventType type;
	Event *next;
};

// Global params
static Graph *graph = NULL;
static IDnum maxtrans = 10;

// Global pointers
static RecycleBin *transcriptMemory = NULL;
static RecycleBin *eventMemory = NULL;

// DEBUG
static Node *heavyNode = NULL;
static boolean debug = false;

Transcript *allocateTranscript()
{
	if (transcriptMemory == NULL)
		transcriptMemory =
		    newRecycleBin(sizeof(Transcript), BLOCK_SIZE);

	return allocatePointer(transcriptMemory);
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

static Coordinate getTotalCoverageOases(Node * node)
{
	Category cat;
	Coordinate res = 0;
	PassageMarkerI marker;

	for (cat = 0; cat < CATEGORIES; cat++)
		res += getVirtualCoverage(node, cat);

	for (marker = getMarker(node); marker;
	     marker = getNextInNode(marker))
		res += getPassageMarkerLength(marker);

	return res;
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
		    && getTotalCoverageOases(locus->contigs[index]) > maxCov) {
			maxCov = getTotalCoverageOases(locus->contigs[index]);
			nextNode = locus->contigs[index];
		}
	}

	if (debug)
		velvetLog("Next node will be %li\n", (long) getNodeID(nextNode));

	heavyNode = nextNode;
	return nextNode;
}

static double bestSuccessorWeight(Node * node, double * scores) {
	Connection *connect;
	double maxWeight = 0;
	double edgeWeight = 0;

	for (connect = getConnection(node); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getTwinNode(getConnectionDestination(connect)))) {
			if (heavyNode && getConnectionDestination(connect) == heavyNode) {
				return 100000000 * getConnectionWeight(connect);
			// DEBUG: modified for real DP
			} else if (scores[getNodeID(getTwinNode(getConnectionDestination(connect))) + nodeCount(graph)] + getConnectionWeight(connect) > maxWeight) {
				maxWeight = scores[getNodeID(getTwinNode(getConnectionDestination(connect))) + nodeCount(graph)] + getConnectionWeight(connect);
				edgeWeight = getConnectionWeight(connect);
			}
		}
	}

	//if (heavyNode && getTwinNode(node) == heavyNode)
	//	edgeWeight *= 100000000;

	return edgeWeight;
}

static Node *chooseBestPredecessor(Node * node, double * scores)
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

		// DEBUG: Modified for real DP
		if (getNodeStatus(getConnectionDestination(connect)) == 2
		    && scores[getNodeID(getConnectionDestination(connect)) + nodeCount(graph)] + getConnectionWeight(connect) > maxWeight) {
			nextNode = getConnectionDestination(connect);
			maxWeight = scores[getNodeID(getConnectionDestination(connect)) + nodeCount(graph)] + getConnectionWeight(connect);
		}
	}

	return nextNode;
}

static boolean detectLoop(Node * start, Node * current, double * scores, IDnum maxlength)
{
	Node * node = chooseBestPredecessor(current, scores);

	if (debug)
	    velvetLog("Testing nodes %li -> %li for loop (%li, %li)\n", (long) getNodeID(current), (long) getNodeID(node), (long) getNodeID(start), (long) maxlength);

	if (node == NULL || maxlength == 0)
		return false;
	else if (node == getTwinNode(start)) {
		if (debug)
		    puts("LOOP");
		return true;
	} else
		return detectLoop(start, node, scores, maxlength - 1);
}

static void computeDPScore(Node * node, double *scores, IDnum contigCount)
{
	Connection *connect;
	double maxWeight = 0;
	double edgeWeight = 0;
	Node *predecessor = NULL;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect))) {
			if (heavyNode && getConnectionDestination(connect) == heavyNode) {
				maxWeight = scores[getNodeID(getConnectionDestination(connect)) + nodeCount(graph)] + getConnectionWeight(connect) * 100000000;
				edgeWeight = getConnectionWeight(connect);
				predecessor = getConnectionDestination(connect);
				break;
			} else if (scores[getNodeID(getConnectionDestination(connect)) + nodeCount(graph)] + getConnectionWeight(connect) > maxWeight) {
				maxWeight = scores[getNodeID(getConnectionDestination(connect)) + nodeCount(graph)] + getConnectionWeight(connect);
				edgeWeight = getConnectionWeight(connect);
				predecessor = getConnectionDestination(connect);
			}
		}
	}

	if (predecessor == NULL) {
	    if (heavyNode && node == heavyNode)
		scores[getNodeID(node) + nodeCount(graph)] = 100000000;
	    else
		scores[getNodeID(node) + nodeCount(graph)] = 0;
	    if (debug)
		velvetLog("Score: %li -> %f\n", (long) getNodeID(node), scores[getNodeID(node) + nodeCount(graph)]);
	    return;
	}

	if (heavyNode && node == heavyNode) 
		scores[getNodeID(node) + nodeCount(graph)] =
		    100000000 * edgeWeight + scores[getNodeID(predecessor) + nodeCount(graph)];
	else if (getNodeStatus(getTwinNode(node)) == 2 && detectLoop(node, node, scores, contigCount)) {
		if (getUniqueness(node))
		    scores[getNodeID(node) + nodeCount(graph)] = 0;
		else
		    scores[getNodeID(node) + nodeCount(graph)] = 
			    scores[getNodeID(predecessor) + nodeCount(graph)] - bestSuccessorWeight(node, scores);
	} else 
		scores[getNodeID(node) + nodeCount(graph)] = maxWeight;

	if (debug)
	    velvetLog("Score: %li -> %f\n", (long) getNodeID(node), scores[getNodeID(node) + nodeCount(graph)]);
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

static Connection *getReverseMarkedConnection(Node * node)
{
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect))
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect)))
			return connect;

	return NULL;
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

static void createAlternativeDPStartFromNode(Node * node, Locus * locus)
{
	Node *current = node;
	Node *next;

	setSingleNodeStatus(current, 3);

	while ((next = chooseBestUnvisitedPredecessor(current))) {
		current = next;
		setSingleNodeStatus(current, 3);
	}

	if (getConnectionBetweenNodes(node, getTwinNode(current))) {
		destroyNodeBackConnections(node);
		simplifyFromNode(node, locus);
	} else {
		destroyNodeBackConnections(current);
		simplifyFromNode(current, locus);
	}
}

static boolean lookForPossibleAlternativeDPStartFromNode(Node * node, Locus * locus)
{
	Connection *connect;

	if (getNodeStatus(node) != 1)
		return false;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect)) ==
		    2) {
			createAlternativeDPStartFromNode(node, locus);
			return true;
		}
	}

	return false;
}

static boolean forcePossibleAlternativeDPStartFromNode(Node * node, Locus * locus)
{
	Connection *connect, *connect2;

	if (getNodeStatus(node) != 1)
		return false;

	for (connect = getConnection(getTwinNode(node)); connect != NULL;
	     connect = getNextConnection(connect)) {
		if (getConnectionStatus(connect)
		    && getNodeStatus(getConnectionDestination(connect))) {
			for (connect2 = getConnection(getTwinNode(node)); connect2 != NULL;
			     connect2 = getNextConnection(connect2)) 
				if (getNodeStatus(getConnectionDestination(connect2)))
					setConnectionStatus(connect2, false);
			simplifyFromNode(node, locus);
			//velvetLog("Forced new start %li\n",(long)  getNodeID(node));
			return true;
		}
	}

	return false;
}

static boolean createAlternativeDPStarts(Locus * locus)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++) {
		if (lookForPossibleAlternativeDPStartFromNode
		    (locus->contigs[index], locus)) {
			renumberLocusNodes(locus);
			return true;
		}
	}

	for (index = 0; index < locus->contigCount; index++) {
		if (forcePossibleAlternativeDPStartFromNode
		    (locus->contigs[index], locus)) {
			renumberLocusNodes(locus);
			return true;
		}
	}


	return false;
}

static IDnum extractMajorityPath(Node * maxNode, double * scores)
{
	Node *node = maxNode;
	IDnum nodesInList = 1;

	recordNode(node);
	setSingleNodeStatus(node, 1);

	while (getReverseMarkedConnection(node)) {
		if (debug)
		    velvetLog("Back tracking through node %li\n", (signed long) getNodeID(node));
		node = chooseBestPredecessor(node, scores);
		if (node == NULL)
			break;
		recordNode(node);
		setSingleNodeStatus(node, 1);
		nodesInList++;
		if (scores[getNodeID(node) + nodeCount(graph)] == 0)
			break;
	}
	if (debug)
	    velvetLog("Arrived in %li\n", (long) getNodeID(node));

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
		//velvetLog("Back tracking through node %li\n", (signed long) getNodeID(node));
		node = getConnectionDestination(connect);
		recordNode(node);
		setNodeStatus(node, false);
		nodesInList++;
	}
	//velvetLog("Arrived in %li\n", (long) getNodeID(node));

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
			    getNodeLength(node)/2;
			transcript->distances[index - 1] -=
			    getNodeLength(transcript->contigs[index - 1])/2;
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
	if (!existsMarkedNode()) {
		createAlternativeDPStartFromNode(locus->contigs[0], locus);
		computeHighestExpressedLocusTranscript(locus, scores);
		return;
	}
	// Propagate DP
	while ((node = popNodeRecord())) {
		//velvetLog("Visiting node %li\n", (long) getNodeID(node));
		setSingleNodeStatus(node, 2);
		nodesToVisit--;

		computeDPScore(node, scores, locus->contigCount);
		if (scores[getNodeID(node) + nodeCount(graph)] > overallMaxScore
		    && (!heavyNode || node != getTwinNode(heavyNode))) {
			overallMaxScore = scores[getNodeID(node) + nodeCount(graph)];
			maxNode = node;
		}
		propagateDP(node);

		// If loop prevents nodes from being visited:
		if (nodesToVisit > 0 && !existsMarkedNode()) {
			while (popNodeRecord());
			if (createAlternativeDPStarts(locus))
				computeHighestExpressedLocusTranscript
				    (locus, scores);
			return;
		}
	}

	nodesInList = extractMajorityPath(maxNode, scores);
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

	//velvetLog("%li nodes explained out of %li\n", (long) counter, (long) locus->contigCount);

	return counter;
}

static void makeTranscriptOfNode(Locus * locus, Node* node) {
	Transcript *transcript = allocateTranscript();
	transcript->contigCount = 1;
	transcript->contigs = callocOrExit(1, Node *);
	transcript->contigs[0] = node;
	transcript->distances = callocOrExit(1, Coordinate);
	//transcript->confidence = 1 / (double) locus->contigCount;
	transcript->confidence = 1;
	transcript->next = locus->transcript;
	locus->transcript = transcript;
}

static void expressUnexplainedLongContigs(Locus * locus) {
	IDnum index;
	Node *node;
	Transcript *transcript;

	// Clear the scene
	setLocusStatus(locus, false);

	// Mark used contigs
	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next) {
		for (index = 0; index < transcript->contigCount; index++) {
			node = transcript->contigs[index];
			setSingleNodeStatus(node, true);
		}
	}

	// Find unused contigs
	for (index = 0; index < locus->contigCount; index++) {
		node = locus->contigs[index];
		if (!getNodeStatus(node))
			makeTranscriptOfNode(locus, node);
	}
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

	if (counter >= maxtrans)
		expressUnexplainedLongContigs(locus);

	setLocusStatus(locus, true);
	setLocusConnectionStatus(locus, false);
	setLocusStatus(locus, false);
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
	destroyRecycleBin(eventMemory);
	transcriptMemory = NULL;
	eventMemory = NULL;
	cleanScaffoldMemory();
	cleanNodeListMemory();
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
						       scores);
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
		locus = &(loci[index]);
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

void detachShortReads(ReadSet * reads, int wordLength)
{
	IDnum index;
	IDnum pairID;
	IDnum sequenceCount = reads->readCount;
	IDnum *mateReads = reads->mateReads;

	if (mateReads == NULL)
		return;

	for (index = 0; index < sequenceCount; index++) {
		if (getLength(getTightStringInArray(reads->tSequences, index)) >= wordLength || reads->categories[index] % 2 == 0 )
			continue;

		if (isSecondInPair(reads, index))
		    pairID = index - 1;
		else
		    pairID = index + 1;

		reads->categories[index] = (reads->categories[index] / 2) * 2;
		reads->categories[pairID] = (reads->categories[pairID] / 2) * 2;
	}
}
