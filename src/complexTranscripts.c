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
#include "transcript.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

// Global params
static Graph *graph = NULL;
static IDnum maxtrans = 10;

// DEBUG
static Node *heavyNode = NULL;
static boolean debug = false;

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

	for (transcript = getTranscript(locus); transcript;
	     transcript = getNextTranscript(transcript))
		for (index = 0; index < getTranscriptContigCount(transcript); index++)
			setSingleNodeStatus(getTranscriptContig(transcript, index),
					    true);

	for (index = 0; index < getContigCount(locus); index++) {
		if (!getNodeStatus(getContig(locus, index))
		    && getTotalCoverageOases(getContig(locus, index)) > maxCov) {
			maxCov = getTotalCoverageOases(getContig(locus, index));
			nextNode = getContig(locus, index);
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

	for (index = 0; index < getContigCount(locus); index++) {
		if (lookForPossibleAlternativeDPStartFromNode
		    (getContig(locus, index), locus)) {
			renumberLocusNodes(locus);
			return true;
		}
	}

	for (index = 0; index < getContigCount(locus); index++) {
		if (forcePossibleAlternativeDPStartFromNode
		    (getContig(locus, index), locus)) {
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

void computeHighestExpressedLocusTranscript(Locus * locus, double *scores)
{
	IDnum index;
	Node *node;
	double overallMaxScore = -1;
	Node *maxNode = NULL;
	IDnum nodesToVisit = getContigCount(locus);
	IDnum nodesInList;

	// Clear node statuses:
	// Record all nodes from which to start the DP
	setLocusStatus(locus, true);

	for (index = 0; index < getContigCount(locus); index++) {
		node = getContig(locus, index);
		if (getReverseMarkedConnection(node) == NULL)
			recordNode(node);
	}

	// OK, let's name a volunteer...
	if (!existsMarkedNode()) {
		createAlternativeDPStartFromNode(getContig(locus, 0), locus);
		computeHighestExpressedLocusTranscript(locus, scores);
		return;
	}
	// Propagate DP
	while ((node = popNodeRecord())) {
		//velvetLog("Visiting node %li\n", (long) getNodeID(node));
		setSingleNodeStatus(node, 2);
		nodesToVisit--;

		computeDPScore(node, scores, getContigCount(locus));
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
	findHeaviestNonUsedNode(locus);
}

static IDnum explainedContigs(Locus * locus)
{
	IDnum index;
	Node *node;
	Transcript *transcript;
	IDnum counter = 0;

	setLocusStatus(locus, false);

	for (transcript = getTranscript(locus); transcript != NULL;
	     transcript = getNextTranscript(transcript)) {
		for (index = 0; index < getTranscriptContigCount(transcript); index++) {
			node = getTranscriptContig(transcript, index);
			if (!getNodeStatus(node)) {
				counter++;
				setSingleNodeStatus(node, true);
			}
		}
	}

	//velvetLog("%li nodes explained out of %li\n", (long) counter, (long) getContigCount(locus));

	return counter;
}

static void makeTranscriptOfNode(Locus * locus, Node* node) {
	Transcript *transcript = newTranscript(1, 1);
	addContigToTranscript(transcript, node, 0);
	addTranscript(locus, transcript);
}

static void expressUnexplainedLongContigs(Locus * locus) {
	IDnum index;
	Node *node;
	Transcript *transcript;

	// Clear the scene
	setLocusStatus(locus, false);

	// Mark used contigs
	for (transcript = getTranscript(locus); transcript != NULL;
	     transcript = getNextTranscript(transcript)) {
		for (index = 0; index < getTranscriptContigCount(transcript); index++) {
			node = getTranscriptContig(transcript, index);
			setSingleNodeStatus(node, true);
		}
	}

	// Find unused contigs
	for (index = 0; index < getContigCount(locus); index++) {
		node = getContig(locus, index);
		if (!getNodeStatus(node))
			makeTranscriptOfNode(locus, node);
	}
}

void computeHighlyExpressedLocusTranscripts(Locus * locus,
						   double *scores, Graph * argGraph)
{
	IDnum counter = 0;
	graph = argGraph;

	heavyNode = NULL;

	setLocusStatus(locus, true);
	setLocusConnectionStatus(locus, true);

	while (counter++ < maxtrans
	       && explainedContigs(locus) < getContigCount(locus))
		computeHighestExpressedLocusTranscript(locus, scores);

	if (counter >= maxtrans)
		expressUnexplainedLongContigs(locus);

	setLocusStatus(locus, true);
	setLocusConnectionStatus(locus, false);
	setLocusStatus(locus, false);
}
