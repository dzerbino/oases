/*
Copyright 2007, 2008 Daniel Zerbino (zerbino@ebi.ac.uk)

    This file is part of Velvet.

    Velvet is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Velvet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Velvet; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

*/
#include <stdlib.h>
#include <stdio.h>

#include "globals.h"
#include "graph.h"
#include "tightString.h"
#include "dfibHeap.h"
#include "recycleBin.h"
#include "passageMarker.h"
#include "transcript.h"
#include "locallyCorrectedGraph.h"
#include "utility.h"
#include "scaffold.h"
#include "transcript.h"

static const Time INDEL = 0;
static const Time SIM[4][4] = {
	{1, 0, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 1, 0},
	{0, 0, 0, 1}
};

//Global variables used throughout this procedure(internal use only !)
static int MAXREADLENGTH = 100;
static int MAXNODELENGTH = 200;
static double MAXDIVERGENCE = 0.2;
static int MAXGAPS = 3;

static Time *times = NULL;
static Node **previous = NULL;

static DFibHeapNode **dheapNodes = NULL;
static DFibHeap *dheap = NULL;

static TightString *fastSequence = NULL;
static TightString *slowSequence = NULL;

static int WORDLENGTH = 0;
static int SELF_LOOP_CUTOFF = 0;
static Graph *graph = NULL;
static Node *start = NULL;

static PassageMarkerI fastPath = NULL_IDX;
static PassageMarkerI slowPath = NULL_IDX;

static double **Fmatrix;
//End of global variables;

static void setNodeTime(Node * node, Time time)
{
	times[getNodeID(node) + nodeCount(graph)] = time;
}

static Time getNodeTime(Node * node)
{
	return times[getNodeID(node) + nodeCount(graph)];
}

static Node *getNodePrevious(Node * node)
{
	return previous[getNodeID(node) + nodeCount(graph)];
}

static boolean isPreviousToNode(Node * previous, Node * target)
{
	Node *currentNode = target;
	Node *previousNode = NULL;
	Time targetTime = getNodeTime(target);

	//printf("Testing if %li is previous to %li\n", getNodeID(previous), getNodeID(target));

	while (true) {
		//printf("CCC %li %f\n", getNodeID(currentNode), getNodeTime(currentNode));

		if (currentNode == previous)
			return true;

		if (currentNode == previousNode)
			return false;

		if (getNodeID(currentNode) > nodeCount(graph)
		    || getNodeID(currentNode) < -nodeCount(graph)) {
			printf("Node ID??? %d %d\n",
			       getNodeID(currentNode),
			       getNodeID(previousNode));
		}

		if (getNodeTime(currentNode) != targetTime)
			return false;

		previousNode = currentNode;
		currentNode = getNodePrevious(currentNode);
	}
}

static boolean
extractSequence(PassageMarkerI path, TightString * sequence)
{
	PassageMarkerI marker;
	Coordinate seqLength = 0;
	Coordinate writeIndex = 0;

	//printf("Extracting sequence %li ... ", pathLength);

	//Measure length
	for (marker = getNextInSequence(path); !isTerminal(marker);
	     marker = getNextInSequence(marker))
		seqLength += getNodeLength(getNode(marker));

	if (seqLength > MAXREADLENGTH)
		return false;
	else
		setTightStringLength(sequence, seqLength);

	//Copy sequences
	for (marker = getNextInSequence(path); !isTerminal(marker);
	     marker = getNextInSequence(marker)) {
		appendNodeSequence(getNode(marker), sequence, writeIndex);
		writeIndex += getNodeLength(getNode(marker));
	}

	return true;
}

static Time max(Time A, Time B, Time C)
{
	if (A >= B && A >= C)
		return A;
	else if (B >= C)
		return B;
	else
		return C;
}

static boolean
compareSequences(TightString * sequence1, TightString * sequence2)
{
	Coordinate i, j;
	Coordinate length1 = getLength(sequence1);
	Coordinate length2 = getLength(sequence2);
	Coordinate maxLength;
	Time Choice1, Choice2, Choice3;
	Time maxScore;

	if (length1 == 0 || length2 == 0)
		return false;

	maxLength = (length1 > length2 ? length1 : length2);

	if (length1 < WORDLENGTH || length2 < WORDLENGTH)
		if (maxLength - length1 > MAXGAPS
		    || maxLength - length2 > MAXGAPS)
			return false;

	for (i = 0; i <= length1; i++)
		Fmatrix[i][0] = 0;
	for (j = 0; j <= length2; j++)
		Fmatrix[0][j] = 0;

	for (i = 1; i <= length1; i++) {
		for (j = 1; j <= length2; j++) {
			Choice1 =
			    Fmatrix[i - 1][j - 1] +
			    SIM[(int) getNucleotide(i - 1, sequence1)]
			    [(int) getNucleotide(j - 1, sequence2)];
			Choice2 = Fmatrix[i - 1][j] + INDEL;
			Choice3 = Fmatrix[i][j - 1] + INDEL;
			Fmatrix[i][j] = max(Choice1, Choice2, Choice3);
		}
	}

	maxScore = Fmatrix[length1][length2];

	if ((1 - maxScore / maxLength) > MAXDIVERGENCE)
		return false;

	return true;
}

static void destroyPaths()
{
	PassageMarkerI marker;

	while (slowPath != NULL_IDX) {
		marker = slowPath;
		getNodeTime(getNode(marker));
		getNodeTime(getTwinNode(getNode(marker)));

		slowPath = getNextInSequence(marker);
		destroyPassageMarker(marker);
	}

	while (fastPath != NULL_IDX) {
		marker = fastPath;
		getNodeTime(getNode(marker));
		getNodeTime(getTwinNode(getNode(marker)));
		fastPath = getNextInSequence(marker);
		destroyPassageMarker(marker);
	}
}

static void cleanUpRedundancy_local2()
{
	PassageMarkerI current;

	for (current = getNextInSequence(slowPath); !isTerminal(current);
	     current = getNextInSequence(current))
		removeNodeFromLocus(getNode(current));

	destroyPaths();
}

static void comparePaths_local2(Node * destination, Node * origin)
{
	IDnum slowLength, fastLength;
	Node *fastNode, *slowNode;
	IDnum i;
	PassageMarkerI marker;

	//Measure lengths
	slowLength = fastLength = 0;
	fastNode = destination;
	slowNode = origin;

	//puts("Looking into separate paths");

	while (fastNode != slowNode) {
		//printf("Fast node %li Slow node %li\n", getNodeID(fastNode), getNodeID(slowNode));

		if (getNodeTime(fastNode) > getNodeTime(slowNode)) {
			fastLength++;
			fastNode = getNodePrevious(fastNode);
		} else if (getNodeTime(fastNode) < getNodeTime(slowNode)) {
			slowLength++;
			slowNode = getNodePrevious(slowNode);
		} else if (isPreviousToNode(slowNode, fastNode)) {
			while (fastNode != slowNode) {
				fastLength++;
				fastNode = getNodePrevious(fastNode);
			}
		} else if (isPreviousToNode(fastNode, slowNode)) {
			while (slowNode != fastNode) {
				slowLength++;
				slowNode = getNodePrevious(slowNode);
			}
		} else {
			fastLength++;
			fastNode = getNodePrevious(fastNode);
			slowLength++;
			slowNode = getNodePrevious(slowNode);
		}

		if (slowLength > MAXNODELENGTH
		    || fastLength > MAXNODELENGTH) {
			//printf("Paths too fragmented %li %li\n", slowLength, fastLength);
			return;
		}
	}

	if (fastLength == 0)
		return;

	//Backtracking to record actual paths
	fastPath = addUncertainPassageMarker(1, destination);
	setPassageMarkerStatus(fastPath, true);

	for (i = 0; i < fastLength; i++) {
		marker =
		    addUncertainPassageMarker(1,
					      getNodePrevious(getNode
							      (fastPath)));
		setPassageMarkerStatus(marker, true);
		connectPassageMarkers(marker, fastPath, graph);
		fastPath = marker;
	}

	slowPath = addUncertainPassageMarker(2, destination);
	setPassageMarkerStatus(slowPath, true);

	marker = addUncertainPassageMarker(2, origin);
	setPassageMarkerStatus(marker, true);
	connectPassageMarkers(marker, slowPath, graph);
	slowPath = marker;

	for (i = 0; i < slowLength; i++) {
		marker =
		    addUncertainPassageMarker(2,
					      getNodePrevious(getNode
							      (slowPath)));
		setPassageMarkerStatus(marker, true);
		connectPassageMarkers(marker, slowPath, graph);
		slowPath = marker;
	}

	//Extract sequences
	if (!extractSequence(fastPath, fastSequence)
	    || !extractSequence(slowPath, slowSequence)) {
		//puts("Paths too long");
		destroyPaths();
		return;
	}
	//Compare sequences
	if (compareSequences(fastSequence, slowSequence)) {
		//puts("Correcting discrepancy");
		cleanUpRedundancy_local2();
		return;
	}
	//puts("\tFinished comparing paths, changes made");
	destroyPaths();
}

static void tourBusConnection_local2(Node * origin, Connection * connect, Time originTime)
{
	Node *destination = getTwinNode(getConnectionDestination(connect));
	Time connectTime, totalTime, destinationTime;
	IDnum nodeIndex = getNodeID(destination) + nodeCount(graph);
	Node *oldPrevious = previous[nodeIndex];

	//printf("Trying connection from %li -> %li\n", getNodeID(origin), getNodeID(destination)); 

	if (oldPrevious == origin)
		return;

	connectTime =
	    ((Time) getNodeLength(origin)) / ((Time) getConnectionWeight(connect));
	totalTime = originTime + connectTime;

	destinationTime = times[nodeIndex];

	if (destinationTime == -1) {
		//puts("New destination");
		setNodeTime(destination, totalTime);
		dheapNodes[nodeIndex] =
		    insertNodeIntoDHeap(dheap, totalTime, destination);
		previous[nodeIndex] = origin;
		return;
	} else if (destinationTime > totalTime) {
		//printf("Previously visited from slower node %li\n", getNodeID(getNodePrevious(destination))); 
		if (dheapNodes[nodeIndex] == NULL) {
			return;
		}

		setNodeTime(destination, totalTime);
		replaceKeyInDHeap(dheap, dheapNodes[nodeIndex], totalTime);
		previous[nodeIndex] = origin;

		comparePaths_local2(destination, oldPrevious);
		return;
	} else {
		//printf("Previously visited by faster node %li\n", getNodeID(getNodePrevious(destination))); 
		comparePaths_local2(destination, origin);
	}
}

static void tourBusNode_local2(Node * node)
{
	Connection * connect;
	Node *destination;
	Time nodeTime = getNodeTime(node);

	//printf("Node %li %f %i %p\n", getNodeID(node),
	//       times[getNodeID(node) + nodeCount(graph)], simpleArcCount(node),
	//       node);

	for (connect = getConnection(node); connect != NULL; connect = getNextConnection(connect)) {
		destination = getTwinNode(getConnectionDestination(connect));

		// Node doesn't belong to the marked node area 
		if (getNodeStatus(destination) != 1)
			continue;

		tourBusConnection_local2(node, connect, nodeTime);

		if (getNodeStatus(node) != 1)
			break;
	}
}

static void tourBus_local2(Node * startingPoint)
{
	Node *currentNode = startingPoint;
	IDnum nodeID = getNodeID(startingPoint) + nodeCount(graph);

	//printf("Tour bus from node %li...\n", getNodeID(startingPoint));

	times[nodeID] = 0;
	previous[nodeID] = currentNode;

	while (currentNode != NULL) {
		dheapNodes[getNodeID(currentNode) + nodeCount(graph)] =
		    NULL;
		tourBusNode_local2(currentNode);
		currentNode = removeNextNodeFromDHeap(dheap);
	}
}

void prepareGraphForLocalCorrections2(Graph * argGraph)
{
	IDnum nodes = nodeCount(argGraph);
	IDnum index;

	//Setting global params
	graph = argGraph;
	WORDLENGTH = getWordLength(graph);;
	SELF_LOOP_CUTOFF = WORDLENGTH;
	// Done with global params

	printf("Preparing to correct graph with cutoff %f\n",
	       MAXDIVERGENCE);

	// Allocating memory
	times = mallocOrExit(2 * nodes + 1, Time);
	previous = mallocOrExit(2 * nodes + 1, Node *);

	dheapNodes = mallocOrExit(2 * nodes + 1, DFibHeapNode *);

	printf("Allocating array %p\n", dheapNodes);

	dheap = newDFibHeap();

	fastSequence = newTightString(MAXREADLENGTH);
	slowSequence = newTightString(MAXREADLENGTH);

	for (index = 0; index < (2 * nodeCount(graph) + 1); index++) {
		times[index] = -1;
		dheapNodes[index] = NULL;
		previous[index] = NULL;
	}

	Fmatrix = callocOrExit(MAXREADLENGTH + 1, double *);
	for (index = 0; index < MAXREADLENGTH + 1; index++)
		Fmatrix[index] = callocOrExit(MAXREADLENGTH + 1, double);
	//Done with memory 
}

void correctGraphLocally2(Node * argStart, Locus * locus)
{
	IDnum index, nodeIndex;

	start = argStart;
	//printf("Correcting graph from node %li\n", getNodeID(start));

	for (index = 0; index < locus->contigCount; index++) {
		nodeIndex = getNodeID(locus->contigs[index]) + nodeCount(graph);
		times[nodeIndex] = -1;
		dheapNodes[nodeIndex] = NULL;
		previous[nodeIndex] = NULL;
	}

	tourBus_local2(start);
}

void deactivateLocalCorrectionSettings2()
{
	puts("Deactivating local correction settings");
	IDnum index;

	for (index = 0; index <= MAXREADLENGTH; index++) {
		free(Fmatrix[index]);
	}
	free(Fmatrix);

	free(times);
	free(previous);
	free(dheapNodes);
	destroyDHeap(dheap);

	destroyTightString(fastSequence);
	destroyTightString(slowSequence);
}

void setLocalMaxReadLength2(int value)
{
	MAXREADLENGTH = value;
	MAXNODELENGTH = 2 * value;
}

void setLocalMaxGaps2(int value)
{
	MAXGAPS = value;
}

void setLocalMaxDivergence2(double value)
{
	MAXDIVERGENCE = value;
}
