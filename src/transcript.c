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
#include "transcript.h"

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
static NodeList *markedNodes;
static RecycleBin *nodeListMemory = NULL;
static RecycleBin *transcriptMemory = NULL;
static RecycleBin *eventMemory = NULL;

// DEBUG
static Node *heavyNode = NULL;
static boolean debug = false;

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

	//velvetLog("Recording node %li\n", (long) getNodeID(node));

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

	//velvetLog("Popping record %li\n", (long) getNodeID(node));

	deallocateNodeList(nodeList);
	return node;
}

Transcript *allocateTranscript()
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
	//velvetLog("Extending from node %li\n", (long) getNodeID(node));

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

#ifndef COLOR
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
#else 
	if (getNodeLength(node) < 2)
		return 0;

	// Initialise:
	for (index = 0; index < 6; index++)
		totalStops[index] = 0;
	n2 = getNucleotideInNode(node, 0);

	// Scna sequence and spot Stop codons
	for (index = 1; index < nodeLength; index++) {
		n1 = n2;
		n2 = getNucleotideInNode(node, index);

		if (n1 == THYMINE 
		    && (n2 == ADENINE || n2 == GUANINE))
			totalStops[index % 3]++;

		if (n1 == CYTOSINE && n2 == GUANINE)
			totalStops[index % 3]++;

		if (n2 == THYMINE 
		    && (n1 == ADENINE || n1 == GUANINE))
			totalStops[3 + index % 3]++;

		if (n1 == GUANINE && n2 == CYTOSINE)
			totalStops[3 + index % 3]++;
	}

#endif

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

static void setNodeConnectionStatus(Node * node, boolean status)
{
	Connection *connect;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			setConnectionStatus(connect, status);
}

void removeNodeFromLocus(Node * node) {
	setNodeStatus(node, false);
	setNodeConnectionStatus(node, false);
}

static void renumberLocusNodes(Locus * locus) {
	IDnum index;
	Node * node;
	IDnum counter = 0;
	Node ** newArray;

	for (index = 0; index < locus->contigCount; index++) {
		node = locus->contigs[index];
		if (!getNodeStatus(node)) {
			locus->contigs[index] = NULL;
			counter++;
			if (getUniqueness(node))
				locus->longContigCount--;	
		}
	}
	
	if (counter == 0)
		return;

	newArray = callocOrExit(locus->contigCount - counter, Node *);
	counter = 0;	

	for (index = 0; index < locus->contigCount; index++) {
		node = locus->contigs[index];

		if (node == NULL) 
			counter++;
		else 
			newArray[index - counter] = node;
	}

	free(locus->contigs);
	locus->contigs = newArray;
	locus->contigCount -= counter;
}

static Connection *getReverseActiveConnection(Node * node)
{
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus(getConnectionDestination(connect)))
			return connect;

	return NULL;
}

static void simplifyFromNode(Node * node, Locus * locus) {
	if (getReverseActiveConnection(node))
		return;

	correctGraphLocally2(node, locus);
}

static void setLocusStatus(Locus * locus, boolean status)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		setSingleNodeStatus(locus->contigs[index], status);
}

static void simplifyLocus(Locus * locus) {
	IDnum index;
	
	setLocusStatus(locus, true);

	for (index = 0; index < locus->contigCount; index++) 
		simplifyFromNode(locus->contigs[index], locus);

	renumberLocusNodes(locus);
	setLocusStatus(locus, false);
}

static void simplifyLoci(Locus * loci, IDnum locusCount) {
	IDnum index;

	prepareGraphForLocalCorrections2(graph);

	for (index = 0; index < locusCount; index++) 
		simplifyLocus(&(loci[index]));

	deactivateLocalCorrectionSettings2();
}

Locus *extractGraphLoci(Graph * argGraph, ReadSet * reads,
			boolean * dubious, ShortLength * lengths,
			IDnum * locusCount, boolean scaffolding)
{
	Locus *loci;

	graph = argGraph;

	buildScaffold(graph, reads, dubious, lengths, scaffolding);

	puts("Extracting loci from connection graph...");

	*locusCount = countConnectedComponents(graph);
	velvetLog("Counted %li mRNA loci\n", (long) *locusCount);

	loci = extractConnectedComponents(*locusCount);
	if (doubleStrandedGraph(graph))
		orientLoci(loci, *locusCount);

	transitiveReduction();
	simplifyLoci(loci, *locusCount);

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
	if (markedNodes == NULL) {
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
		if (nodesToVisit > 0 && markedNodes == NULL) {
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
	destroyRecycleBin(nodeListMemory);
	destroyRecycleBin(eventMemory);
	nodeListMemory = NULL;
	transcriptMemory = NULL;
	eventMemory = NULL;
	cleanScaffoldMemory();
}

static Coordinate getTranscriptLength(Transcript * transcript) {
	IDnum index;
	Coordinate totalLength = 0;

	if (transcript->contigCount == 0)
		return 0;

	totalLength = getWordLength(graph) - 1;

	for (index = 0; index < transcript->contigCount; index++) {
		totalLength += getNodeLength(transcript->contigs[index]);
		if (index < transcript->contigCount - 1) {
			if (transcript->distances[index] < getWordLength(graph) && getNodeLength(transcript->contigs[index+1]) >= getWordLength(graph) - 1) {
				totalLength += transcript->distances[index];
			} else if (transcript->distances[index] > 0 && transcript->distances[index] < 10) {
				totalLength += 10;
			} else {
				totalLength += transcript->distances[index];
			}
		}
	}

	return totalLength;
}

static char revComp(char base) {
	if (base == 'A')
		return 'T';
	if (base == 'T')
		return 'A';
	if (base == 'G')
		return 'C';
	if (base == 'C')
		return 'G';
	return 'N';
}

static char * nSequence(Coordinate length) {
	char * sequence = callocOrExit(length + 1, char);
	Coordinate position;

	for (position = 0; position < length; position++) 
		sequence[position] = 'N';
	sequence[length] = '\0';

	return sequence;
}

static void printFastA(FILE * outfile, char * sequence) {
	Coordinate position;
	size_t length = strlen(sequence);
	char str[100];

	for (position = 0; position < length; position += strlen(str)) {
		strncpy(str, sequence + position, 60);
		str[60] = '\0';
		fprintf(outfile, "%s\n", str);
	}
}

static void addLongNodeToSequence(char * sequence, Node * node, Coordinate offset) {
	Coordinate position;
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	int wordShift = getWordLength(graph) - 1;

	for (position = 0; position < strlen(string); position++) 
	    if (sequence[offset + position - wordShift] == 'N')
		sequence[offset + position - wordShift] = string[position]; 
	free(string);
}

static void addShortNodeToSequence(char * sequence, Node * node, Coordinate offset) {
	Coordinate position;
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	int wordShift = getWordLength(graph) - 1;

	for (position = 0; position < strlen(string); position++) 
	    if (sequence[offset + position] == 'N')
		sequence[offset + position] = string[position]; 

	free(string);
	string = expandNodeFragment(getTwinNode(node), 0, getNodeLength(node),
				      getWordLength(graph));

	for (position = 0; position < strlen(string); position++) 
	    if (sequence[offset - wordShift + getNodeLength(node) - 1 - position] == 'N')
		sequence[offset - wordShift + getNodeLength(node) - 1 - position] = revComp(string[position]); 
	free(string);
}

static void addNodeToSequence(char * sequence, Node * node, Coordinate offset) {
	int wordShift = getWordLength(graph) - 1;

	// Add sequence
	if (getNodeLength(node) >= wordShift) 
		addLongNodeToSequence(sequence, node, offset);
	else
		addShortNodeToSequence(sequence, node, offset);
}

static void exportTranscript(Transcript * transcript, IDnum locusID,
			     IDnum transID, IDnum transcriptCount, FILE * outfile)
{
	IDnum index;
	Coordinate offset = getWordLength(graph) - 1;
	char * sequence = nSequence(getTranscriptLength(transcript)); 

	// Extract sequence
	for (index = 0; index < transcript->contigCount; index++) {
	    addNodeToSequence(sequence, transcript->contigs[index], offset);

	    // Increment offset
	    offset += getNodeLength(transcript->contigs[index]);
	    if (index < transcript->contigCount - 1)
		    offset += transcript->distances[index];
	}

	// Count initial N's
	for (offset = 0; offset < strlen(sequence); offset++)
		if (sequence[offset] != 'N')
			break;
	
	// Print
	fprintf(outfile, ">Locus_%li_Transcript_%li/%li_Confidence_%.3f_Length_%li\n",
		(long) locusID + 1, (long) transID + 1, (long) transcriptCount, transcript->confidence, (long) (strlen(sequence) - offset));
	printFastA(outfile, sequence + offset);

	free(sequence);
}

static void exportLocusTranscripts(Locus * locus, IDnum locusID,
				   FILE * outfile, Coordinate minTransLength)
{
	IDnum index = 0;
	Transcript *transcript;
	IDnum transcriptCount = 0;

	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next)
		if (getTranscriptLength(transcript) >= minTransLength)
			transcriptCount++;

	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next)
		if (getTranscriptLength(transcript) >= minTransLength)
			exportTranscript(transcript, locusID, index++, transcriptCount, outfile);
}

void exportTranscripts(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength, Graph * argGraph)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;
	graph = argGraph;

	velvetLog("Exporting transcripts to %s\n", filename);

	for (index = 0; index < locusCount; index++)
		exportLocusTranscripts(&(loci[index]), index, outfile, minTransLength);
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

void computeTranscripts(Locus * loci, IDnum locusCount) {
	IDnum index;
	Locus *locus;
	IDnum distribution[4];
	int configuration;
	double *scores = callocOrExit(2 * nodeCount(graph), double);

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
	fprintf(outfile, ">Locus_%li_Node_%li\n", (long) index + 1,
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

	fprintf(outfile, "Locus %li: ", (long) index + 1);

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
		velvetLog("Exporting AS events to %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	for (index = 0; index < locusCount; index++)
		exportLocusASEvents(index, &(loci[index]), outfile);

	fclose(outfile);
}

static void exportTranscriptContigs(Transcript * transcript, IDnum locusID,
				    IDnum transID, IDnum transcriptCount, FILE * outfile)
{
	IDnum index;
	Coordinate totalLength = 0;

	if (transcript->contigCount == 0)
		return;

	// Header
	fprintf(outfile, ">Locus_%li_Transcript_%li/%li_Confidence_%.3f_Length_%li\n",
		(long) locusID + 1, (long) transID + 1, (long) transcriptCount, transcript->confidence, (long) getTranscriptLength(transcript));

	totalLength = getWordLength(graph) - 1;

	// Sequence
	for (index = 0; index < transcript->contigCount; index++) {
		totalLength += getNodeLength(transcript->contigs[index]);
		fprintf(outfile, "%li:%lli",
			(long) getNodeID(transcript->contigs[index]),
			(long long) totalLength);
		if (index < transcript->contigCount - 1) {
			if (transcript->distances[index] < getWordLength(graph) && getNodeLength(transcript->contigs[index+1]) >= getWordLength(graph) - 1) {
				fprintf(outfile, "-(0)->");
				totalLength += transcript->distances[index];
			} else if (transcript->distances[index] > 0 && transcript->distances[index] < 10) {
				fprintf(outfile, "-(10)->");
				totalLength += 10;
			} else {
				fprintf(outfile, "-(%li)->",
					(long) transcript->distances[index]);
				totalLength += transcript->distances[index];
			}
		}
	}

	fprintf(outfile, "\n");
}

static void exportLocusContigs(IDnum locusID, Locus * locus,
			       FILE * outfile, Coordinate minTransLength)
{
	IDnum index = 0;
	Transcript *transcript;
	IDnum transcriptCount = 0;

	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next)
		if (getTranscriptLength(transcript) >= minTransLength)
			transcriptCount++;

	exportLocusNodes(locusID, locus, outfile);
	for (transcript = locus->transcript; transcript != NULL;
	     transcript = transcript->next)
		if (getTranscriptLength(transcript) >= minTransLength)
			exportTranscriptContigs(transcript, locusID, index++, transcriptCount,
						outfile);
}

void exportContigOrders(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;

	if (outfile)
		velvetLog("Exporting transcript contigs to %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	for (index = 0; index < locusCount; index++)
		exportLocusContigs(index, &(loci[index]), outfile, minTransLength);

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
				fprintf(file, "%lli -> %lli [label=%lli_%li_%li",
					-nodeID,
					(long long)
					-getNodeID(getConnectionDestination
						   (connect)),
					(long long)
					getConnectionDistance(connect),
					(long) getConnectionPairedCount(connect),
					(long) getConnectionDirectCount(connect));
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
		velvetLog("Reading read set file %s;\n", filename);
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
	velvetLog("%d sequences found\n", sequenceCount);

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

static void exportAMOSLib(FILE * outfile, Graph * graph, Category cat)
{
	Coordinate distance = getInsertLength(graph, cat * 2);
	double variance = getInsertLength_var(graph, cat * 2);

	if (distance == -1)
		return;

	fprintf(outfile, "{LIB\n");
	fprintf(outfile, "iid:%d\n", (int) (cat + 1));
	fprintf(outfile, "{DST\n");
	fprintf(outfile, "mea:%lld\n", (long long) distance);
	fprintf(outfile, "std:%lld\n", (long long) sqrt(variance));
	fprintf(outfile, "}\n");
	fprintf(outfile, "}\n");
}

static void exportAMOSMarker(FILE * outfile, PassageMarkerI marker,
			     Coordinate nodeLength, Coordinate offset,
			     int wordShift, boolean firstUniqueMet)
{
	Coordinate sequenceStart, sequenceFinish;

	sequenceStart = getPassageMarkerStart(marker);
	sequenceFinish = getPassageMarkerFinish(marker);

	if (firstUniqueMet == -1) {
		if (getPassageMarkerSequenceID(marker) < 0) {
			sequenceFinish += wordShift;
			sequenceStart += wordShift;
		}
	} else if (firstUniqueMet == 0) {
		if (getPassageMarkerSequenceID(marker) > 0) 
			sequenceFinish += wordShift;
		else 
			sequenceStart += wordShift;
	} else if (firstUniqueMet == 1) {
		if (getPassageMarkerSequenceID(marker) > 0) {
			sequenceFinish += wordShift;
			sequenceStart += wordShift;
		}
	}

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%d\n", getAbsolutePassMarkerSeqID(marker));
	fprintf(outfile, "off:%lld\n", (long long) (offset + getStartOffset(marker)));
	fprintf(outfile, "clr:%lld,%lld\n", (long long) sequenceStart, (long long) sequenceFinish);
	fprintf(outfile, "}\n");
}

static void exportAMOSShortMarker(FILE * outfile, ShortReadMarker * marker,
				  ReadSet * reads, Coordinate offset, boolean firstUniqueMet, int wordShift)
{
	Coordinate read_offset =
	    getShortReadMarkerPosition(marker) -
	    getShortReadMarkerOffset(marker)
	    + offset;
	TightString *sequence = getTightStringInArray(reads->tSequences, getShortReadMarkerID(marker) - 1);

	if (firstUniqueMet == 1)
		read_offset -= wordShift;

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%d\n", getShortReadMarkerID(marker));
	fprintf(outfile, "off:%lld\n", (long long) read_offset);
	fprintf(outfile, "clr:0,%lld\n", (long long) getLength(sequence));
	fprintf(outfile, "}\n");
}

static void exportAMOSReverseShortMarker(FILE * outfile,
					 ShortReadMarker * marker,
					 Coordinate nodeLength,
					 int wordShift, ReadSet * reads,
					 Coordinate offset, boolean firstUniqueMet)
{
	TightString *sequence = getTightStringInArray(reads->tSequences, getShortReadMarkerID(marker) - 1);

	Coordinate read_offset =
	    nodeLength - getShortReadMarkerPosition(marker) +
	    getShortReadMarkerOffset(marker) - getLength(sequence) +
	    wordShift + offset;

	if (firstUniqueMet == 1)
		read_offset -= wordShift;

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%d\n", getShortReadMarkerID(marker));
	fprintf(outfile, "off:%lld\n", (long long) read_offset);
	fprintf(outfile, "clr:%lld,0\n", (long long) getLength(sequence));
	fprintf(outfile, "}\n");
}


static void exportAMOSContig(FILE * outfile, ReadSet * reads, Transcript * transcript, 
			     IDnum startIndex, IDnum finishIndex, Graph * graph, 
			     IDnum locusID, IDnum transcriptID, IDnum internalIndex, IDnum iid)
{
	Coordinate start;
	PassageMarkerI marker;
	ShortReadMarker *shortMarkerArray, *shortMarker;
	Coordinate index, maxIndex;
	int wordShift = getWordLength(graph) - 1;
	Coordinate offset = 0;
	IDnum nodeIndex;
	Node * node;
	boolean firstUniqueMet = false;
	Coordinate length = 0;
	int column = 0;
	Nucleotide nucleotide;

	fprintf(outfile, "{CTG\n");
	fprintf(outfile, "iid:%d\n", iid);
	fprintf(outfile, "eid:%d-%d-%d\n", locusID, transcriptID, internalIndex);

	fprintf(outfile, "seq:\n");
	for (nodeIndex = startIndex; nodeIndex <= finishIndex; nodeIndex++) {
		node = transcript->contigs[nodeIndex];
		
		// If first unique met, print initial k-mer
		if (!firstUniqueMet && getNodeLength(node) >= wordShift) {
			for (start = 0; start < wordShift; start++) {
				nucleotide = getNucleotideInNode(getTwinNode(node), getNodeLength(node) - start - 1);
#ifndef COLOR
				nucleotide = 3 - nucleotide;
#endif
				switch (nucleotide) {
				case ADENINE:
					fprintf(outfile, "A");
					break;	
				case CYTOSINE:
					fprintf(outfile, "C");
					break;	
				case GUANINE:
					fprintf(outfile, "G");
					break;	
				case THYMINE:
					fprintf(outfile, "T");
					break;	
				default:
					abort();
				}

				if (column++ == 60) {
					fprintf(outfile, "\n");
					column = 0;
				}
			}
			
			length += wordShift;
			firstUniqueMet = true;
		}

		// Print proper sequence
		if (!firstUniqueMet) {
			for (start = 0; start < getNodeLength(node); start++) {
				nucleotide = getNucleotideInNode(getTwinNode(node), getNodeLength(node) - start - 1);
#ifndef COLOR
				nucleotide = 3 - nucleotide;
#endif
				switch (nucleotide) {
				case ADENINE:
					fprintf(outfile, "A");
					break;	
				case CYTOSINE:
					fprintf(outfile, "C");
					break;	
				case GUANINE:
					fprintf(outfile, "G");
					break;	
				case THYMINE:
					fprintf(outfile, "T");
					break;	
				}

				if (column++ == 60) {
					fprintf(outfile, "\n");
					column = 0;
				}
			}
		} else {
			for (start = 0; start < getNodeLength(node); start++) {
				nucleotide = getNucleotideInNode(node, start);

				switch (nucleotide) {
				case ADENINE:
					fprintf(outfile, "A");
					break;	
				case CYTOSINE:
					fprintf(outfile, "C");
					break;	
				case GUANINE:
					fprintf(outfile, "G");
					break;	
				case THYMINE:
					fprintf(outfile, "T");
					break;	
				}

				if (column++ == 60) {
					fprintf(outfile, "\n");
					column = 0;
				}
			}
		}

		length += getNodeLength(node);
	}
	fprintf(outfile, "\n.\n");

	fprintf(outfile, "qlt:\n");
	column = 0;
	for (start = 0; start < length; start++) {
		fprintf(outfile, "I");

		if (column++ == 60) {
			fprintf(outfile, "\n");
			column = 0;
		}
	}
	fprintf(outfile, "\n.\n");

	firstUniqueMet = -1;
	for (nodeIndex = startIndex; nodeIndex <= finishIndex; nodeIndex++) {
		node = transcript->contigs[nodeIndex];
		if (firstUniqueMet == 0)
			firstUniqueMet = 1;
		else if (firstUniqueMet == -1 && getNodeLength(node) >= wordShift) 
			firstUniqueMet = 0;

		for (marker = getMarker(node); marker != NULL_IDX;
		     marker = getNextInNode(marker))
			exportAMOSMarker(outfile, marker, getNodeLength(node),
					 offset, wordShift, firstUniqueMet);
		

		if (readStartsAreActivated(graph)) {
			shortMarkerArray = getNodeReads(node, graph);
			maxIndex = getNodeReadCount(node, graph);
			for (index = 0; index < maxIndex; index++) {
				shortMarker =
				    getShortReadMarkerAtIndex(shortMarkerArray,
							      index);
				exportAMOSShortMarker(outfile, shortMarker, reads,
						      offset, firstUniqueMet, wordShift);
			}

			shortMarkerArray = getNodeReads(getTwinNode(node), graph);
			maxIndex = getNodeReadCount(getTwinNode(node), graph);
			for (index = 0; index < maxIndex; index++) {
				shortMarker =
				    getShortReadMarkerAtIndex(shortMarkerArray,
							      index);
				exportAMOSReverseShortMarker(outfile, shortMarker,
							     getNodeLength(node),
							     wordShift, reads,
							     offset, firstUniqueMet);
			}
		}

		offset += getNodeLength(node);
		if (firstUniqueMet == 0) 
			offset += wordShift;
	}

	fprintf(outfile, "}\n");
}

static void exportAMOSTranscript(FILE* outfile, ReadSet * reads, Transcript * transcript, IDnum locusID, IDnum transcriptID, Graph * graph) 
{
	IDnum smallIndex = 0;
	static IDnum iid = 1;
	IDnum contigIndex = iid;
	int wordShift = getWordLength(graph) - 1;
	IDnum nodeIndex;
	boolean uniqueBunch = false;
	IDnum startIndex = 0;
	Coordinate start;
	Coordinate contigLength;
	Node * node;

	contigLength = 0;
	for (nodeIndex = 0; nodeIndex < transcript->contigCount; nodeIndex++) {
		node = transcript->contigs[nodeIndex];
		contigLength += getNodeLength(node);
		if (getNodeLength(node) >= wordShift)
			uniqueBunch = true;

		if (nodeIndex == transcript->contigCount - 1 || transcript->distances[nodeIndex] > 0) {
			if (contigLength > 0) {
				exportAMOSContig(outfile, reads, transcript, startIndex, nodeIndex, graph, 
							 locusID, transcriptID, smallIndex++, iid++);
				startIndex = nodeIndex + 1;
			}
			uniqueBunch = false;
			contigLength = 0;
		} 
	}

	fprintf(outfile, "{SCF\n");
	fprintf(outfile, "eid:%d-%d\n", locusID, transcriptID);

	uniqueBunch = false;
	start = 0;
	contigLength = 0;
	for (nodeIndex = 0; nodeIndex < transcript->contigCount; nodeIndex++) {
		node = transcript->contigs[nodeIndex];
		contigLength += getNodeLength(node);
		if (getNodeLength(node) >= wordShift)
			uniqueBunch = true;

		if (nodeIndex == transcript->contigCount - 1 || transcript->distances[nodeIndex] > 0) {
			if (contigLength > 0) {
				if (uniqueBunch)
					contigLength += wordShift;

				fprintf(outfile, "{TLE\n");
				fprintf(outfile, "off:%lld\n", (long long) start);
				fprintf(outfile, "clr:0,%lld\n", (long long) contigLength);
				fprintf(outfile, "src:%d\n", contigIndex++);
				fprintf(outfile, "}\n");
			}
			start += contigLength;
			start += transcript->distances[nodeIndex];
			uniqueBunch = false;
			contigLength = 0;
		} 
	}

	fprintf(outfile, "}\n");
}

static void exportAMOSLocus(FILE * outfile, ReadSet * reads, Locus * locus, IDnum locusID, Coordinate minTransLength,
			   Graph * graph)
{
	Transcript * transcript;
	IDnum transcriptID = 0;

	for (transcript = locus->transcript; transcript; transcript = transcript->next) 
		if (getTranscriptLength(transcript) >= minTransLength)
			exportAMOSTranscript(outfile, reads, transcript, locusID, transcriptID++, graph);
}

static void exportAMOSRead(FILE * outfile, TightString * tString,
			   IDnum index, IDnum frg_index)
{
	Coordinate start, finish;
	char str[100];

	fprintf(outfile, "{RED\n");
	fprintf(outfile, "iid:%d\n", index);
	fprintf(outfile, "eid:%d\n", index);
	if (frg_index > 0)
		fprintf(outfile, "frg:%d\n", frg_index);

	fprintf(outfile, "seq:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "qlt:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "}\n");
}

void exportAMOSTranscripts(Graph * graph,
		       Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char * directory)
{
	IDnum index;
	Category cat;
	Locus * locus;
	FILE *outfile;
	char filename[10000];

	strcpy(filename, directory);
	strcat(filename, "/oases_asm.afg");
	velvetLog("Writing into AMOS file %s...\n", filename);
	outfile = fopen(filename, "w");

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to AMOS file %s",
		       filename);

	for (cat = 0; cat <= CATEGORIES; cat++)
		exportAMOSLib(outfile, graph, cat);

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			fprintf(outfile, "{FRG\n");
			fprintf(outfile, "lib:%d\n",
				(int) ((reads->categories[index - 1] / 2) + 1));
			fprintf(outfile, "rds:%d,%d\n", index,
				index + 1);
			fprintf(outfile, "eid:%d\n", index);
			fprintf(outfile, "iid:%d\n", index);
			fprintf(outfile, "typ:I\n");
			fprintf(outfile, "}\n");
			index++;
		}
	}

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       index);
			index++;
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       index - 1);
		} else {
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       -1);
		}
	}

	for (index = 0; index < locusCount; index++) {
		locus = &(loci[index]);

		if (locus == NULL)
			continue;

		exportAMOSLocus(outfile, reads, locus, index, minTransLength, graph);
	}

	fclose(outfile);

}

static void markUsedReads(Node * node, boolean * used)
{
	IDnum readID;
	ShortReadMarker * shortReadArray, * shortReadMarker;
	IDnum shortReadCount, shortReadIndex;
	PassageMarkerI marker;

	if (node == NULL || getNodeStatus(node))
		return;
	else
		setNodeStatus(node, true);
	
	// Long reads
	for(marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker)) {
		readID = getPassageMarkerSequenceID(marker);
		if (readID < 0)
			readID = -readID;
		used[readID] = true;	
	}	

	// Short reads		
	if (!readStartsAreActivated(graph))
		return;

	shortReadArray = getNodeReads(node, graph);
	shortReadCount = getNodeReadCount(node, graph);
	for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
		shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
		readID = getShortReadMarkerID(shortReadMarker);
		used[readID] = true;	
	}
	
	shortReadArray = getNodeReads(getTwinNode(node), graph);
	shortReadCount = getNodeReadCount(getTwinNode(node), graph);
	for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
		shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
		readID = getShortReadMarkerID(shortReadMarker);
		used[readID] = true;	
	}

}

void exportUnusedTranscriptReads(Graph* graph, Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char* directory) {
	char *outFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	FILE * outfile;
	boolean * used = callocOrExit(sequenceCount(graph) + 1, boolean);
	IDnum nodeID, readID;
	IDnum locusIndex;
	Transcript * transcript;

	strcpy(outFilename, directory);
	strcat(outFilename, "/UnusedReads.fa");
	outfile = fopen(outFilename, "w");

	velvetLog("Printing unused reads into %s\n", outFilename);

	resetNodeStatus(graph);

	for (locusIndex = 0; locusIndex < locusCount; locusIndex++)
		for (transcript = loci[locusIndex].transcript; transcript; transcript = transcript->next)
			if (getTranscriptLength(transcript) >= minTransLength)
				for(nodeID = 0; nodeID < transcript->contigCount; nodeID++)
					markUsedReads(transcript->contigs[nodeID], used);

	for (readID = 1; readID <= sequenceCount(graph); readID++) 
		if (!used[readID])
			exportTightString(outfile, getTightStringInArray(reads->tSequences, readID - 1), readID);	

	free(outFilename);
	free(used);	
	fclose(outfile);
}

IDnum usedTranscriptReads(Graph * graph, Coordinate minTransLength, Locus * loci, IDnum locusCount) 
{
	boolean * used = callocOrExit(sequenceCount(graph) + 1, boolean);
	IDnum nodeID, readID;
	IDnum locusIndex;
	Transcript * transcript;
	IDnum res = 0;

	resetNodeStatus(graph);

	for (locusIndex = 0; locusIndex < locusCount; locusIndex++)
		for (transcript = loci[locusIndex].transcript; transcript; transcript = transcript->next)
			if (getTranscriptLength(transcript) >= minTransLength)
				for(nodeID = 0; nodeID < transcript->contigCount; nodeID++)
					markUsedReads(transcript->contigs[nodeID], used);

	for (readID = 1; readID <= sequenceCount(graph); readID++) 
		if (used[readID])
			res++;

	free(used);	

	return res;
}

void setPairedThreshold(double pairedThreshold) {
	scaffold_setPairedThreshold(pairedThreshold);
}

void setDegreeCutoff(int val) {
	scaffold_setDegreeCutoff(val);
}

static IDnum getReferenceCount(ReadSet * reads) {
	IDnum index;

	for (index = 0; index < reads->readCount; index++) 
		if (reads->categories[index] <= 2 * CATEGORIES + 1)
			break;

	return index;
}

//////////////////////////////////////////////////////////////////////////
// Reference identifiers
//////////////////////////////////////////////////////////////////////////

typedef struct referenceCoord_st ReferenceCoord;

struct referenceCoord_st {
	char * name;
	Coordinate start;
	Coordinate finish;
	boolean positive_strand;
};

static ReferenceCoord * collectReferenceCoords(char * sequencesFilename, IDnum referenceCount) {
	FILE * file = fopen(sequencesFilename, "r");
	char line[MAXLINE];
	char name[500];
	Coordinate start, finish;
	long long longlongvar;
	IDnum refIndex = 0;
	ReferenceCoord * refCoords = callocOrExit(referenceCount, ReferenceCoord);
	int i;

	while (fgets(line, MAXLINE, file)) {
		if (line[0] == '>') {
			if (strchr(line, ':')) {
				sscanf(strtok(line, ":-\r\n"), ">%s", name);
				sscanf(strtok(NULL, ":-\r\n"), "%lli", &longlongvar);
				start = longlongvar;
				sscanf(strtok(NULL, ":-\r\n"), "%lli", &longlongvar);
				finish = longlongvar;
				refCoords[refIndex].name = callocOrExit(strlen(name) + 1, char);  
				if (start <= finish) {
					strcpy(refCoords[refIndex].name, name);
					refCoords[refIndex].start = start;
					refCoords[refIndex].finish = finish;
					refCoords[refIndex].positive_strand = true;
				} else {
					strcpy(refCoords[refIndex].name, name);
					refCoords[refIndex].start = finish;
					refCoords[refIndex].finish = start;
					refCoords[refIndex].positive_strand = false;
				}
			} else {
				for (i = strlen(line) - 1;
				     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
					line[i] = '\0';
				}

				refCoords[refIndex].name = callocOrExit(strlen(name) + 1, char);  
				strcpy(name, line + 1);
				strcpy(refCoords[refIndex].name, name);
				refCoords[refIndex].start = 1;
				refCoords[refIndex].finish = -1;
				refCoords[refIndex].positive_strand = true;
			}
			if (++refIndex == referenceCount)
				break;	
		}
	}
	
	fclose(file);
	return refCoords;
}

typedef struct refMap_st {
	Coordinate start;
	Coordinate finish;
	IDnum refID;
	Coordinate refStart;
	Coordinate refFinish;
} ReferenceMapping; 

static int compareReferenceMappings(const void * A, const void * B) {
	ReferenceMapping * refMapA = (ReferenceMapping *) A;
	ReferenceMapping * refMapB = (ReferenceMapping *) B;
	
	if (refMapA->start < refMapB->start)
		return -1;
	else if (refMapA->start == refMapB->start)
		return 0;
	else 
		return 1;
}

static void initializeReferenceMapping(ReferenceMapping * refMap, PassageMarkerI marker, Transcript * transcript, IDnum nodeIndex, Coordinate nodeOffset) {
	PassageMarkerI finishMarker = marker;
	Coordinate totalLength = getNodeLength(transcript->contigs[nodeIndex]);
	IDnum index;

	for (index = nodeIndex + 1; index < transcript->contigCount; index++) {
		if (!getNextInSequence(finishMarker))
			break;

		if (getNode(getNextInSequence(finishMarker)) == transcript->contigs[index]) {
			finishMarker = getNextInSequence(finishMarker);
			totalLength += getNodeLength(transcript->contigs[index]);
		}
	}

	refMap->start = nodeOffset + getStartOffset(marker);
	refMap->finish = nodeOffset + totalLength - getFinishOffset(finishMarker); 
	refMap->refID = getPassageMarkerSequenceID(marker);
	refMap->refStart = getPassageMarkerStart(marker);
	refMap->refFinish = getPassageMarkerFinish(finishMarker);
}

static void fprintfReferenceMapping(FILE * file, ReferenceMapping * mapping, ReferenceCoord * refCoords, int wordLength) {
	ReferenceCoord * refCoord;
	Coordinate start, finish;

	if (mapping->refID > 0) 
		refCoord = &refCoords[mapping->refID - 1];
	else
		refCoord = &refCoords[-mapping->refID - 1];

	if (mapping->refID > 0) {
		if (refCoord->positive_strand) {
			start = refCoord->start + mapping->refStart;
			finish = refCoord->start + mapping->refFinish + wordLength - 2;
		} else {
			start = refCoord->finish - mapping->refStart + wordLength - 1;
			finish = refCoord->finish - mapping->refFinish + 1;
		}
	} else {
		if (refCoord->positive_strand) {
			start = refCoord->start + mapping->refStart + wordLength - 1;
			finish = refCoord->start + mapping->refFinish + 1;
		} else {
			start = refCoord->finish - mapping->refStart; 
			finish = refCoord->finish - mapping->refFinish + wordLength;  
		}
	}
		
	fprintf(file, "%lli\t%lli\t%s\t%lli\t%lli\n",
		(long long) mapping->start + 1, (long long) mapping->finish + wordLength - 1, 
		refCoord->name, (long long) start, (long long) finish);
}

static void exportTranscriptMapping(FILE * outfile, Transcript * transcript, IDnum locusIndex, IDnum transcriptIndex, ReadSet * reads, ReferenceCoord * refCoords, int wordLength) {
	PassageMarkerI marker;
	ReferenceMapping * referenceMappings;
	IDnum index;
	IDnum referenceCount = 0;
	IDnum nodeIndex;
	Coordinate nodeOffset = 0;

	// Count reference sequences
	for (nodeIndex = 0; nodeIndex < transcript->contigCount; nodeIndex++)
		for (marker = getMarker(transcript->contigs[nodeIndex]); marker; marker = getNextInNode(marker))
			if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] > 2 * CATEGORIES + 1
			    && (nodeIndex == 0 
				|| getNode(getPreviousInSequence(marker)) != transcript->contigs[nodeIndex - 1]))
				referenceCount++;

	// Header
	fprintf(outfile, ">Locus_%li_Transcript_%li\n", (long) locusIndex + 1, (long) transcriptIndex + 1);

	// Create table
	referenceMappings = callocOrExit(referenceCount, ReferenceMapping);	

	// Initialize table
	referenceCount = 0;
	for (nodeIndex = 0; nodeIndex < transcript->contigCount; nodeIndex++) {
		for (marker = getMarker(transcript->contigs[nodeIndex]); marker; marker = getNextInNode(marker))
			if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] > 2 * CATEGORIES + 1
			    && (nodeIndex == 0 
				|| getNode(getPreviousInSequence(marker)) != transcript->contigs[nodeIndex - 1]))
				initializeReferenceMapping(&referenceMappings[referenceCount++], marker, transcript, nodeIndex, nodeOffset);
		nodeOffset += getNodeLength(transcript->contigs[nodeIndex]);
		nodeOffset += transcript->distances[nodeIndex];
	}

	// Sort table
	qsort(referenceMappings, referenceCount, sizeof(ReferenceMapping), compareReferenceMappings);

	// Print table
	for (index = 0; index < referenceCount; index++)
		fprintfReferenceMapping(outfile, &referenceMappings[index], refCoords, wordLength);

	// Clean table
	free(referenceMappings);
}

static void exportLocusMapping(FILE * outfile, Locus * loci, IDnum locusIndex, ReadSet * reads, ReferenceCoord * refCoords, Coordinate minTransLength, int wordLength) {
	Transcript * transcript;
	IDnum transcriptIndex = 0;

	for (transcript = loci[locusIndex].transcript; transcript != NULL;
	     transcript = transcript->next)
		if (getTranscriptLength(transcript) >= minTransLength)
			exportTranscriptMapping(outfile, transcript, locusIndex, transcriptIndex++, reads, refCoords, wordLength);
		
}

void exportTranscriptMappings(Locus * loci, IDnum locusCount, 
			      Graph * graph, ReadSet * reads,
			      Coordinate minLength, char * directory)
{
	FILE * outfile;
	IDnum locusIndex, refIndex;
	ReferenceCoord * refCoords;
	IDnum referenceCount = getReferenceCount(reads); 
	char filename[10000];

	strcpy(filename, directory);
	strcat(filename, "/transcript-alignments.psa");
	velvetLog("Writing into pseudo-alignment file %s...\n", filename);
	outfile = fopen(filename, "w");

	if (referenceCount == 0)	
		return;

	strcpy(filename, directory);
	strcat(filename, "/Sequences");
	refCoords = collectReferenceCoords(filename, referenceCount);

	strcpy(filename, directory);
	strcat(filename, "/transcript-alignments.psa");
	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		velvetLog("Could not write into %s, sorry\n", filename);
		return;
	} else {
		velvetLog("Writing contigs into %s...\n", filename);
	}

	for (locusIndex = 0; locusIndex < locusCount; locusIndex++) 
		exportLocusMapping(outfile, loci, locusIndex, reads, refCoords, minLength, getWordLength(graph));

	for (refIndex = 0; refIndex < referenceCount; refIndex++)
		free(refCoords[refIndex].name);
	free(refCoords);
	fclose(outfile);
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
