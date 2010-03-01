/*
    Copyright 2009,2010 Daniel Zerbino (dzerbino@soe.ucsc.edu)

    This file is part of Oases.

    Oases is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
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
#include <math.h>

#include "globals.h"
#include "recycleBin.h"
#include "utility.h"
#include "graph.h"
#include "passageMarker.h"
#include "readSet.h"
#include "locallyCorrectedGraph.h"
#include "transcript.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#define BLOCK_SIZE  100000
#define LN2 1.4
#define LENGTHCUTOFF 50

typedef struct connection_st Connection;
typedef struct readOccurence_st ReadOccurence;

struct connection_st {
	Node *destination;
	Connection *next;
	Connection *previous;
	Connection *twin;
	double distance;
	double variance;
	double weight;
	IDnum direct_count;
	IDnum paired_count;
	boolean status;
};

struct readOccurence_st {
	Coordinate position;
	Coordinate offset;
	IDnum nodeID;
};

// Global params
static Graph *graph = NULL;
static IDnum UNRELIABLE_CONNECTION_CUTOFF = 4;

// Global pointers
static Connection **scaffold = NULL;
static RecycleBin *connectionMemory = NULL;

static Connection *allocateConnection()
{
	if (connectionMemory == NULL)
		connectionMemory =
		    newRecycleBin(sizeof(Connection), BLOCK_SIZE);

	return allocatePointer(connectionMemory);
}

static void deallocateConnection(Connection * connect)
{
	deallocatePointer(connectionMemory, connect);
}

boolean getConnectionStatus(Connection * connect)
{
	return connect->status;
}

void setConnectionStatus(Connection * connect, boolean status)
{
	connect->status = status;
	if (connect->twin)
		connect->twin->status = status;
}

Coordinate getConnectionDistance(Connection * connect)
{
	return (Coordinate) connect->distance;
}

static double norm(double X)
{
	return 0.4 * exp(-X * X / 2);
}

static double normInt(double X, double Y)
{
	return (erf(0.7 * Y) - erf(0.7 * X)) / 2;
}

static IDnum expectedNumberOfConnections(IDnum IDA, Connection * connect,
					 IDnum ** counts, Category cat)
{
	Node *A = getNodeInGraph(graph, IDA);
	Node *B = connect->destination;
	IDnum IDB = getNodeID(B);
	double left, middle, right;
	Coordinate longLength, shortLength, D;
	double M, N, O, P;
	Coordinate mu = getInsertLength(graph, cat);
	double sigma = sqrt(getInsertLength_var(graph, cat));
	double result;
	double densityA, densityB, minDensity;

	if (mu <= 0)
		return 0;

	if (getNodeLength(A) == 0 || getNodeLength(B) == 0)
		return 0;

	if (getNodeLength(A) < getNodeLength(B)) {
		longLength = getNodeLength(B);
		shortLength = getNodeLength(A);
	} else {
		longLength = getNodeLength(A);
		shortLength = getNodeLength(B);
	}

	densityA = counts[cat][IDA + nodeCount(graph)] / getNodeLength(A);
	densityB = counts[cat][IDB + nodeCount(graph)] / getNodeLength(B);
	minDensity = densityA > densityB ? densityB : densityA;

	D = getConnectionDistance(connect) - (longLength +
					      shortLength) / 2;

	M = (D - mu) / sigma;
	N = (D + shortLength - mu) / sigma;
	O = (D + longLength - mu) / sigma;
	P = (D + shortLength + longLength - mu) / sigma;

	left = ((norm(M) - norm(N)) - M * normInt(M, N)) * sigma;
	middle = shortLength * normInt(N, O);
	right = ((norm(O) - norm(P)) - P * normInt(O, P)) * (-sigma);

	result = (minDensity * (left + middle + right));

	if (result > 0)
		return (IDnum) result;
	else
		return 0;
}

Connection *getConnection(Node * node)
{
	return scaffold[getNodeID(node) + nodeCount(graph)];
}

void destroyConnection(Connection * connect, IDnum nodeID)
{
	Connection *previous, *next;

	//printf("Destroying connection from %li to %li\n", nodeID, getNodeID(connect->destination));

	if (connect == NULL)
		return;

	previous = connect->previous;
	next = connect->next;

	if (previous != NULL)
		previous->next = next;
	if (next != NULL)
		next->previous = previous;

	if (scaffold[nodeID + nodeCount(graph)] == connect)
		scaffold[nodeID + nodeCount(graph)] = next;

	if (connect->twin != NULL) {
		connect->twin->twin = NULL;
		destroyConnection(connect->twin,
				  getNodeID(connect->destination));
	}

	deallocateConnection(connect);
}

static boolean testConnection(IDnum IDA, Connection * connect,
			      IDnum ** counts)
{
	IDnum total = 0;
	Category cat;

	// Destroy tenuous connections
	if (connect->weight < 0.1)
		return false;

	if (connect->paired_count + connect->direct_count <
	    UNRELIABLE_CONNECTION_CUTOFF)
		return false;

	if (!getUniqueness(connect->destination)
	    && connect->direct_count == 0)
		return false;
	else if (!getUniqueness(connect->destination))
		return true;

	for (cat = 0; cat <= CATEGORIES; cat++)
		total +=
		    expectedNumberOfConnections(IDA, connect, counts, cat);

	if (total == 0 && connect->direct_count == false)
		return false;

	// Remove inconsistent connections
	return connect->paired_count >= total / 10;
}

static IDnum *computeReadToNodeCounts()
{
	IDnum readIndex, nodeIndex;
	IDnum maxNodeIndex = 2 * nodeCount(graph) + 1;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	IDnum *readNodeCounts = callocOrExit(maxReadIndex, IDnum);
	boolean *readMarker = callocOrExit(maxReadIndex, boolean);
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	Node *node;
	IDnum nodeReadCount;

	//puts("Computing read to node mapping array sizes");

	for (nodeIndex = 0; nodeIndex < maxNodeIndex; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodeCount(graph));
		if (node == NULL)
			continue;
		nodeArray = getNodeReads(node, graph);
		nodeReadCount = getNodeReadCount(node, graph);

		// Short reads
		for (readIndex = 0; readIndex < nodeReadCount; readIndex++) {
			shortMarker =
			    getShortReadMarkerAtIndex(nodeArray,
						      readIndex);
			readNodeCounts[getShortReadMarkerID
				       (shortMarker)]++;
		}

		// Long reads
		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex < 0)
				continue;

			if (readMarker[readIndex])
				continue;

			readNodeCounts[readIndex]++;
			readMarker[readIndex] = true;
		}

		// Clean up marker array
		for (marker = getMarker(node); marker != NULL;
		     marker = getNextInNode(marker)) {
			readIndex = getPassageMarkerSequenceID(marker);
			if (readIndex > 0)
				readMarker[readIndex] = false;
		}
	}

	free(readMarker);
	return readNodeCounts;
}

static ReadOccurence **allocateReadToNodeTables(IDnum * readNodeCounts)
{
	IDnum readIndex;
	IDnum maxReadIndex = sequenceCount(graph) + 1;
	ReadOccurence **readNodes =
	    callocOrExit(maxReadIndex, ReadOccurence *);

	for (readIndex = 1; readIndex < maxReadIndex; readIndex++) {
		if (readNodeCounts[readIndex] != 0) {
			readNodes[readIndex] =
			    callocOrExit(readNodeCounts[readIndex],
					 ReadOccurence);
			readNodeCounts[readIndex] = 0;
		}
	}

	return readNodes;
}

static void computePartialReadToNodeMapping(IDnum nodeID,
					    ReadOccurence ** readNodes,
					    IDnum * readNodeCounts,
					    boolean * readMarker)
{
	ShortReadMarker *shortMarker;
	IDnum index, readIndex;
	ReadOccurence *readArray, *readOccurence;
	Node *node = getNodeInGraph(graph, nodeID);
	ShortReadMarker *nodeArray = getNodeReads(node, graph);
	IDnum nodeReadCount = getNodeReadCount(node, graph);
	PassageMarker *marker;

	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		readIndex = getShortReadMarkerID(shortMarker);
		readArray = readNodes[readIndex];
		readOccurence = &readArray[readNodeCounts[readIndex]];
		readOccurence->nodeID = nodeID;
		readOccurence->position =
		    getShortReadMarkerPosition(shortMarker);
		readOccurence->offset =
		    getShortReadMarkerOffset(shortMarker);
		readNodeCounts[readIndex]++;
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex < 0)
			continue;

		if (!readMarker[readIndex]) {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex]];
			readOccurence->nodeID = nodeID;
			readOccurence->position = getStartOffset(marker);
			readOccurence->offset =
			    getPassageMarkerStart(marker);
			readNodeCounts[readIndex]++;
			readMarker[readIndex] = true;
		} else {
			readArray = readNodes[readIndex];
			readOccurence =
			    &readArray[readNodeCounts[readIndex] - 1];
			readOccurence->position = -1;
			readOccurence->offset = -1;
		}
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		readIndex = getPassageMarkerSequenceID(marker);
		if (readIndex > 0)
			readMarker[readIndex] = false;
	}
}

static ReadOccurence **computeReadToNodeMappings(IDnum * readNodeCounts)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	ReadOccurence **readNodes =
	    allocateReadToNodeTables(readNodeCounts);
	boolean *readMarker =
	    callocOrExit(sequenceCount(graph) + 1, boolean);

	//puts("Computing read to node mappings");

	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
		if (nodeID != 0 && getNodeInGraph(graph, nodeID))
			computePartialReadToNodeMapping(nodeID, readNodes,
							readNodeCounts,
							readMarker);

	free(readMarker);
	return readNodes;
}

static Connection *findConnection(IDnum nodeID, IDnum node2ID)
{
	Node *node2 = getNodeInGraph(graph, node2ID);
	Connection *connect;

	if (node2 == NULL)
		return NULL;

	for (connect = scaffold[nodeID + nodeCount(graph)];
	     connect != NULL; connect = connect->next)
		if (connect->destination == node2)
			break;

	return connect;
}

static void createTwinConnection(IDnum nodeID, IDnum node2ID,
				 Connection * connect)
{
	Connection *newConnection = allocateConnection();
	IDnum nodeIndex = nodeID + nodeCount(graph);

	// Fill in
	newConnection->distance = connect->distance;
	newConnection->variance = connect->variance;
	newConnection->direct_count = connect->direct_count;
	newConnection->paired_count = connect->paired_count;
	newConnection->destination = getNodeInGraph(graph, node2ID);
	newConnection->weight = 0;
	newConnection->status = false;

	// Batch to twin
	newConnection->twin = connect;
	connect->twin = newConnection;

	// Insert in scaffold
	newConnection->previous = NULL;
	newConnection->next = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->previous = newConnection;
	scaffold[nodeIndex] = newConnection;
}

static Connection *createNewConnection(IDnum nodeID, IDnum node2ID,
				       IDnum direct_count,
				       IDnum paired_count,
				       Coordinate distance,
				       double variance)
{
	Node *destination = getNodeInGraph(graph, node2ID);
	IDnum nodeIndex = nodeID + nodeCount(graph);
	Connection *connect = allocateConnection();

	// Fill in 
	connect->destination = destination;
	connect->direct_count = direct_count;
	connect->paired_count = paired_count;
	connect->distance = (double) distance;
	connect->variance = variance;
	connect->weight = 0;
	connect->status = false;

	// Insert in scaffold
	connect->previous = NULL;
	connect->next = scaffold[nodeIndex];
	if (scaffold[nodeIndex] != NULL)
		scaffold[nodeIndex]->previous = connect;
	scaffold[nodeIndex] = connect;

	if (node2ID != nodeID)
		createTwinConnection(node2ID, nodeID, connect);
	else 
		connect->twin = NULL;

	return connect;
}

static void readjustConnection(Connection * connect, Coordinate distance,
			       double variance, IDnum direct_count,
			       IDnum paired_count)
{
	connect->direct_count += direct_count;
	connect->paired_count += paired_count;
	connect->distance =
	    (variance * connect->distance +
	     distance * connect->variance) / (variance +
					      connect->variance);
	connect->variance =
	    (variance *
	     connect->variance) / (variance + connect->variance);

	if (connect->twin != NULL) {
		connect->twin->distance = connect->distance;
		connect->twin->variance = connect->variance;
		connect->twin->direct_count = connect->direct_count;
		connect->twin->paired_count = connect->paired_count;
	}
}

static void createConnection(IDnum nodeID, IDnum node2ID,
			     IDnum direct_count,
			     IDnum paired_count,
			     Coordinate distance, double variance)
{
	Connection *connect = findConnection(nodeID, node2ID);

	//printf("Creating connection %li -> %li\n", (long) nodeID, (long) node2ID);

	if (connect != NULL)
		readjustConnection(connect, distance, variance,
				   direct_count, paired_count);
	else
		createNewConnection(nodeID, node2ID, direct_count,
				    paired_count, distance, variance);
}

Connection *getConnectionBetweenNodes(Node * nodeA, Node * nodeB)
{
	Connection *connect;

	for (connect = getConnection(nodeA); connect != NULL;
	     connect = connect->next)
		if (connect->destination == nodeB)
			return connect;

	return NULL;
}

static void incrementConnectionWeight(Connection * connect,
				      double increment)
{
	connect->weight += increment;
	if (connect->weight > 1000)
		connect->weight = 1000;
	if (connect->twin) 
		connect->twin->weight = connect->weight;
}

static void projectFromSingleRead(Node * node,
				  ReadOccurence * readOccurence,
				  Coordinate position,
				  Coordinate offset, Coordinate length,
				  boolean weight)
{
	Coordinate distance = 0;
	Connection *connect;
	Node *target = getNodeInGraph(graph, -readOccurence->nodeID);
	double variance = 1;

	// Filter out troublemakers
	if (readOccurence->position == -1 && readOccurence->offset == -1)
		return;

	if (offset < 0 || readOccurence->offset < 0)
		return;

	if (target == getTwinNode(node) || target == node)
		return;

	if (weight) {
		if ((connect = getConnectionBetweenNodes(node, target)))
			incrementConnectionWeight(connect, 1);
		return;
	}

	if (position < 0) {
		variance += getNodeLength(node) * getNodeLength(node) / 16;
		distance += getNodeLength(node) / 2;
	} else {
		// variance += 0;
		distance += position - offset - getNodeLength(node) / 2;
	}

	if (readOccurence->position < 0) {
		variance +=
		    getNodeLength(target) * getNodeLength(target) / 16;
		distance += getNodeLength(target) / 2;
	} else {
		// variance += 0;
		distance +=
		    -readOccurence->position + readOccurence->offset +
		    getNodeLength(target) / 2;
	}

	if (offset < readOccurence->offset) {
		createConnection(getNodeID(node), getNodeID(target), 1, 0,
				 distance, variance);
	} else {
		createConnection(-getNodeID(node), -getNodeID(target), 1,
				 0, -distance, variance);
	}
}

#define K 0.39894

static void projectFromReadPair(Node * node, ReadOccurence * readOccurence,
				Coordinate position, Coordinate offset,
				Coordinate insertLength,
				double insertVariance, boolean weight)
{
	Coordinate distance = insertLength;
	Coordinate variance = insertVariance;
	Node *target = getNodeInGraph(graph, readOccurence->nodeID);
	Connection *connect;
	double score;

	// Filter for useless reads:
	if (readOccurence->position == -1 && readOccurence->offset == -1)
		return;

	if (target == getTwinNode(node) || target == node)
		return;

	if (getUniqueness(target) && getNodeID(target) < getNodeID(node))
		return;

	if (weight) {
		if (position > 0 && readOccurence->position > 0
		    && (connect =
			getConnectionBetweenNodes(node, target))) {
			distance = getConnectionDistance(connect);
			distance -=
			    position - offset - getNodeLength(node) / 2;
			distance -=
			    readOccurence->position -
			    readOccurence->offset -
			    getNodeLength(target) / 2;
			score =
			    K / sqrt(insertVariance) *
			    exp((insertLength - distance) * (distance -
							     insertLength)
				/ (2 * insertVariance));

			incrementConnectionWeight(connect, score);
		}
		return;
	}

	if (position < 0) {
		variance += getNodeLength(node) * getNodeLength(node) / 16;
		// distance += 0;
	} else {
		// variance += 0;
		distance += position - offset - getNodeLength(node) / 2;
	}

	if (readOccurence->position < 0) {
		variance +=
		    getNodeLength(target) * getNodeLength(target) / 16;
		//distance += 0;
	} else {
		// variance += 0;
		distance +=
		    readOccurence->position - readOccurence->offset -
		    getNodeLength(target) / 2;
	}

	if (distance - getNodeLength(node) / 2 -
	    getNodeLength(target) / 2 < -6 * sqrt(insertVariance))
		return;

	createConnection(getNodeID(node), getNodeID(target), 0, 1,
			 distance, variance);
}

static void projectFromShortRead(Node * node,
				 ShortReadMarker * shortMarker,
				 IDnum * readPairs, Category * cats,
				 ReadOccurence ** readNodes,
				 IDnum * readNodeCounts,
				 Coordinate * lengths, boolean weight)
{
	IDnum index;
	IDnum readIndex = getShortReadMarkerID(shortMarker);
	ReadOccurence *readArray;
	IDnum readPairIndex;
	Category cat;
	Coordinate position = getShortReadMarkerPosition(shortMarker);
	Coordinate offset = getShortReadMarkerOffset(shortMarker);
	Coordinate length = lengths[getShortReadMarkerID(shortMarker) - 1];
	Coordinate insertLength;
	double insertVariance;

	// Filter to remove useless reads
	if (position == -1 && offset == -1)
		return;

	// Going through single-read information
	if (readNodeCounts[readIndex] > 1) {
		readArray = readNodes[readIndex];
		for (index = 0; index < readNodeCounts[readIndex]; index++)
			projectFromSingleRead(node, &readArray[index],
					      position, offset, length,
					      weight);
	}
	// Going through paired read information
	if (readPairs == NULL)
		return;

	readPairIndex = readPairs[readIndex - 1] + 1;

	if (readPairIndex == 0)
		return;

	cat = cats[readIndex - 1];
	insertLength = getInsertLength(graph, cat);
	insertVariance = getInsertLength_var(graph, cat);

	readArray = readNodes[readPairIndex];
	for (index = 0; index < readNodeCounts[readPairIndex]; index++)
		projectFromReadPair(node, &readArray[index], position,
				    offset, insertLength, insertVariance,
				    weight);

}

static void projectFromLongRead(Node * node, PassageMarker * marker,
				IDnum * readPairs, Category * cats,
				ReadOccurence ** readNodes,
				IDnum * readNodeCounts,
				Coordinate * lengths, boolean weight)
{
	IDnum index;
	IDnum readIndex = getPassageMarkerSequenceID(marker);
	ReadOccurence *readArray;
	IDnum readPairIndex;
	Category cat;
	Coordinate position = getStartOffset(marker);
	Coordinate offset = getPassageMarkerStart(marker);
	Coordinate length =
	    lengths[getPassageMarkerSequenceID(marker) - 1];
	Coordinate insertLength;
	double insertVariance;

	// Going through single-read information
	if (readNodeCounts[readIndex] > 1 && position > 0) {
		readArray = readNodes[readIndex];
		for (index = 0; index < readNodeCounts[readIndex]; index++)
			projectFromSingleRead(node, &readArray[index],
					      position, offset, length,
					      weight);
	}
	// Going through paired read information
	if (readPairs == NULL)
		return;

	readPairIndex = readPairs[readIndex - 1] + 1;

	if (readPairIndex == 0)
		return;

	cat = cats[readIndex - 1];
	insertLength = getInsertLength(graph, cat);
	insertVariance = getInsertLength_var(graph, cat);

	readArray = readNodes[readPairIndex];
	for (index = 0; index < readNodeCounts[readPairIndex]; index++)
		projectFromReadPair(node, &readArray[index], position,
				    offset, insertLength, insertVariance,
				    weight);

}

static void projectFromNode(IDnum nodeID,
			    ReadOccurence ** readNodes,
			    IDnum * readNodeCounts,
			    IDnum * readPairs, Category * cats,
			    boolean * dubious,
			    Coordinate * lengths, boolean weight)
{
	IDnum index;
	ShortReadMarker *nodeArray, *shortMarker;
	PassageMarker *marker;
	Node *node;
	IDnum nodeReadCount;

	node = getNodeInGraph(graph, nodeID);

	if (node == NULL || !getUniqueness(node))
		return;

	nodeArray = getNodeReads(node, graph);
	nodeReadCount = getNodeReadCount(node, graph);
	for (index = 0; index < nodeReadCount; index++) {
		shortMarker = getShortReadMarkerAtIndex(nodeArray, index);
		if (dubious[getShortReadMarkerID(shortMarker) - 1])
			continue;
		projectFromShortRead(node, shortMarker, readPairs, cats,
				     readNodes, readNodeCounts, lengths,
				     weight);
	}

	for (marker = getMarker(node); marker != NULL;
	     marker = getNextInNode(marker)) {
		if (getPassageMarkerSequenceID(marker) > 0)
			projectFromLongRead(node, marker, readPairs, cats,
					    readNodes, readNodeCounts,
					    lengths, weight);
	}
}

static Connection **computeNodeToNodeMappings(ReadOccurence ** readNodes,
					      IDnum * readNodeCounts,
					      IDnum * readPairs,
					      Category * cats,
					      boolean * dubious,
					      Coordinate * lengths,
					      boolean weight)
{
	IDnum nodeID;
	IDnum nodes = nodeCount(graph);
	if (!weight)
		scaffold = callocOrExit(2 * nodes + 1, Connection *);

	//puts("Computing direct node to node mappings");

	for (nodeID = -nodes; nodeID <= nodes; nodeID++)
		projectFromNode(nodeID, readNodes, readNodeCounts,
				readPairs, cats, dubious, lengths, weight);

	if (weight)
		return NULL;
	else
		return scaffold;
}

Node *getConnectionDestination(Connection * connect)
{
	return connect->destination;
}

Connection *getNextConnection(Connection * connect)
{
	return connect->next;
}

static void computeLocalNodeToNodeMappingsFromConnections(Connection *
							  connect,
							  Connection *
							  connect2)
{
	Node *node1 = getTwinNode(getConnectionDestination(connect));
	Node *node2 = getTwinNode(getConnectionDestination(connect2));
	IDnum nodeID1 = getNodeID(node1);
	IDnum nodeID2 = getNodeID(node2);
	Coordinate distance =
	    (getNodeLength(node1) + getNodeLength(node2)) / 2;
	Arc *arc;

	if (getUniqueness(node1) || getUniqueness(node2))
		return;


	if ((arc = getArcBetweenNodes(node1, node2, graph))
	    && !getConnectionBetweenNodes(node1, getTwinNode(node2))) {
		createConnection(nodeID1, -nodeID2, getMultiplicity(arc),
				 0, distance,
				 1 / (double) getMultiplicity(arc));
		incrementConnectionWeight(getConnectionBetweenNodes
					  (node1, getTwinNode(node2)),
					  getMultiplicity(arc));
	}

	if ((arc = getArcBetweenNodes(node2, node1, graph))
	    && !getConnectionBetweenNodes(node2, getTwinNode(node1))) {
		createConnection(nodeID2, -nodeID1, getMultiplicity(arc),
				 0, distance,
				 1 / (double) getMultiplicity(arc));
		incrementConnectionWeight(getConnectionBetweenNodes
					  (node2, getTwinNode(node1)),
					  getMultiplicity(arc));
	}

}

static void computeLocalNodeToNodeMappingsFromNode(Node * node)
{
	Connection *connect, *connect2;

	for (connect = getConnection(node); connect;
	     connect = connect->next)
		for (connect2 = getConnection(node); connect2 != connect;
		     connect2 = connect2->next)
			computeLocalNodeToNodeMappingsFromConnections
			    (connect, connect2);
}

static void computeLocalNodeToNodeMappings()
{
	IDnum index;
	Node *node;

	puts("Computing local connections");
	activateArcLookupTable(graph);

	for (index = -nodeCount(graph); index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node && getUniqueness(node))
			computeLocalNodeToNodeMappingsFromNode(node);
	}

	deactivateArcLookupTable(graph);
}

static IDnum **countShortReads(Graph * graph, Category * categories)
{
	IDnum **counts = callocOrExit(CATEGORIES + 1, IDnum *);
	Category cat;
	IDnum nodeIndex;
	IDnum nodes = nodeCount(graph);
	Node *node;
	ShortReadMarker *array, *marker;
	IDnum readCount, readIndex, readID;

	// Allocate memory where needed
	for (cat = 0; cat <= CATEGORIES; cat++)
		if (getInsertLength(graph, cat) > 0)
			counts[cat] =
			    callocOrExit(2 * nodeCount(graph) + 1, IDnum);

	// Start fillin'
	for (nodeIndex = 0; nodeIndex < 2 * nodes + 1; nodeIndex++) {
		node = getNodeInGraph(graph, nodeIndex - nodes);

		if (node == NULL)
			continue;

		array = getNodeReads(node, graph);
		readCount = getNodeReadCount(node, graph);
		for (readIndex = 0; readIndex < readCount; readIndex++) {
			marker =
			    getShortReadMarkerAtIndex(array, readIndex);
			readID = getShortReadMarkerID(marker);
			if (categories)
				cat = categories[readID - 1];
			else
				cat = 1;
			if (cat % 2 == 1 && counts[cat / 2] != NULL)
				counts[cat / 2][nodeIndex]++;
		}
	}

	return counts;
}

static void removeUnreliableConnections(Category * categories)
{
	IDnum maxNodeIndex = nodeCount(graph) * 2 + 1;
	IDnum index;
	Connection *connect, *next;
	Category cat;
	IDnum **counts = countShortReads(graph, categories);
	IDnum nodes = nodeCount(graph);

	for (index = 0; index < maxNodeIndex; index++) {
		for (connect = scaffold[index]; connect != NULL;
		     connect = next) {
			next = connect->next;
			if (!testConnection
			    (index - nodes, connect, counts))
				destroyConnection(connect, index - nodes);
		}
	}

	// Free memory
	for (cat = 0; cat <= CATEGORIES; cat++)
		if (counts[cat])
			free(counts[cat]);
	free(counts);
}

static void defineUniqueness(Node * node)
{
	setUniqueness(node, getNodeLength(node) > LENGTHCUTOFF);
}

static void defineUniqueNodes()
{
	IDnum index;

	for (index = 1; index <= nodeCount(graph); index++)
		defineUniqueness(getNodeInGraph(graph, index));
}

void printOasesConnections(Category * categories)
{
	IDnum maxNodeIndex = nodeCount(graph) * 2 + 1;
	IDnum index;
	Connection *connect, *next;
	Node *node;
	IDnum **counts = countShortReads(graph, categories);
	IDnum nodes = nodeCount(graph);
	Category cat;

	puts("CONNECT IDA IDB dcount pcount dist lengthA lengthB var countA countB coordA coordB real exp distance test");

	for (index = 0; index < maxNodeIndex; index++) {
		node = getNodeInGraph(graph, index - nodeCount(graph));
		for (connect = scaffold[index]; connect != NULL;
		     connect = next) {
			next = connect->next;
			printf
			    ("CONNECT %ld %ld %ld %ld %lld %lld %lld %f %ld %ld",
			     (long) index - nodeCount(graph),
			     (long) getNodeID(connect->destination),
			     (long) connect->direct_count,
			     (long) connect->paired_count,
			     (long long) getConnectionDistance(connect),
			     (long long) getNodeLength(node), (long long)
			     getNodeLength(connect->destination),
			     connect->variance,
			     (long) getNodeReadCount(node, graph),
			     (long) getNodeReadCount(connect->destination,
						     graph));
			if (markerCount(node) == 1
			    && markerCount(connect->destination) == 1)
				printf(" %lld %lld %lld", (long long)
				       getPassageMarkerFinish(getMarker
							      (node)),
				       (long long)
				       getPassageMarkerFinish(getMarker
							      (connect->
							       destination)),
				       (long
					long) (getPassageMarkerFinish
					       (getMarker(node)) -
					       getPassageMarkerFinish
					       (getMarker
						(connect->destination))));
			else
				printf(" ? ? ?");
			printf(" %ld",
			       (long) expectedNumberOfConnections(index -
								  nodeCount
								  (graph),
								  connect,
								  counts,
								  0));
			printf(" %lld",
			       (long long) (getConnectionDistance(connect)
					    - (getNodeLength(node) +
					       getNodeLength
					       (connect->destination)) /
					    2));
			if (testConnection(index - nodes, connect, counts))
				puts(" OK");
			else
				puts(" NG");
		}
	}

	for (cat = 0; cat <= CATEGORIES; cat++)
		if (counts[cat])
			free(counts[cat]);
	free(counts);
}

void setConnectionWeight(Connection * connect, double weight)
{
	connect->weight = weight;
}

double getConnectionWeight(Connection * connect)
{
	return connect->weight;
}

IDnum getConnectionDirectCount(Connection * connect)
{
	return connect->direct_count;
}

// Merges two lists of connections in order of increasing position (used in mergeSort mainly)
static Connection *mergeConnectionLists(Connection * left,
					Connection * right)
{
	Connection *mergedList = NULL;
	Connection *tail = NULL;

	// Choose first element:
	if (left->distance <= right->distance) {
		mergedList = left;
		tail = left;
		left = left->next;
	} else {
		mergedList = right;
		tail = right;
		right = right->next;
	}

	// Iterate while both lists are still non empty
	while (left != NULL && right != NULL) {
		if (left->distance <= right->distance) {
			tail->next = left;
			left->previous = tail;
			left = left->next;
		} else {
			tail->next = right;
			right->previous = tail;
			right = right->next;
		}

		tail = tail->next;
	}

	// Concatenate the remaining list at the end of the merged list
	if (left != NULL) {
		tail->next = left;
		left->previous = tail;
	}

	if (right != NULL) {
		tail->next = right;
		right->previous = tail;
	}

	return mergedList;
}

static void connectionMergeSort(Connection ** connectPtr, IDnum count)
{

	IDnum half = count / 2;
	Connection *left = *connectPtr;
	Connection *ptr = left;
	Connection *right;
	IDnum index;

	if (count == 0 || count == 1)
		return;

	if (count == 2) {
		if ((*connectPtr)->distance >
		    (*connectPtr)->next->distance) {
			(*connectPtr)->next->next = *connectPtr;
			(*connectPtr)->previous = (*connectPtr)->next;
			*connectPtr = (*connectPtr)->next;
			(*connectPtr)->next->next = NULL;
			(*connectPtr)->previous = NULL;
		}
		return;
	}

	for (index = 0; index < half - 1; index++) {
		ptr = ptr->next;
		if (ptr == NULL)
			return;
	}

	right = ptr->next;
	ptr->next = NULL;
	right->previous = NULL;

	connectionMergeSort(&left, half);
	connectionMergeSort(&right, count - half);
	*connectPtr = mergeConnectionLists(left, right);
}

static void sortNodeConnections(IDnum index)
{
	Connection *connect;
	IDnum count = 0;

	for (connect = scaffold[index]; connect != NULL;
	     connect = connect->next)
		count++;

	if (count == 0)
		return;

	connect = scaffold[index];
	connectionMergeSort(&connect, count);

	scaffold[index] = connect;
}

static void sortScaffold()
{
	IDnum index;

	for (index = 0; index <= 2 * nodeCount(graph); index++)
		sortNodeConnections(index);
}

static IDnum countConnections(Node * node)
{
	Connection *connect;
	IDnum count = 0;

	for (connect = getConnection(node); connect;
	     connect = connect->next)
		count++;

	return count;
}

#define VACANT 0
#define INPLAY 1
#define ELIMINATED 2

#define FUZZ 50

void transitiveReductionAtNode(Node * node)
{
	Connection *c, *c2, *c3, *c4;
	Node *dest, *dest2;
	IDnum nodeID = getNodeID(node);
	IDnum connectionCount = countConnections(node);

	for (c = getConnection(node); c; c = c->next)
		setSingleNodeStatus(getTwinNode(c->destination), INPLAY);

	for (c = getConnection(node); c; c = c->next) {
		dest = getTwinNode(c->destination);

		for (c2 = getConnection(dest); c2; c2 = c2->next) {
			dest2 = getTwinNode(c2->destination);

			if (getNodeStatus(dest2) == INPLAY) {
				c3 = getConnectionBetweenNodes(node,
							       c2->
							       destination);

				if (c3->distance >
				    (c->distance + c2->distance) * 0.9) {
					// To avoid circular eliminations
					if ((c4 =
					     getConnectionBetweenNodes
					     (dest2, getTwinNode(dest)))
					    && c->distance >
					    (c3->distance +
					     c4->distance) * 0.9) {
						if (c3->distance <
						    c->distance)
							continue;

						else if (c3->distance ==
							 c->distance
							 &&
							 getNodeStatus
							 (dest) != INPLAY)
							continue;
					}
					setSingleNodeStatus(dest2,
							    ELIMINATED);
					if (--connectionCount <= 1)
						break;
				}
			}
		}

		if (connectionCount <= 1)
			break;
	}

	for (c = getConnection(node); c; c = c2) {
		c2 = c->next;
		dest = getTwinNode(c->destination);
		if (getNodeStatus(dest) == ELIMINATED)
			destroyConnection(c, nodeID);
		setSingleNodeStatus(dest, VACANT);
	}
}

void transitiveReduction()
{
	IDnum index;
	Node *node;

	puts("Transitive reduction of graph");

	sortScaffold();
	resetNodeStatus(graph);

	for (index = -nodeCount(graph); index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (node && getUniqueness(node))
			transitiveReductionAtNode(node);
	}
}

void buildScaffold(Graph * argGraph, ReadSet * reads, boolean * dubious,
		   Coordinate * lengths)
{
	IDnum *readPairs = reads->mateReads;
	Category *cats = reads->categories;
	IDnum *readNodeCounts;
	ReadOccurence **readNodes;
	IDnum index;

	graph = argGraph;

	defineUniqueNodes();
	readNodeCounts = computeReadToNodeCounts();
	readNodes = computeReadToNodeMappings(readNodeCounts);
	scaffold =
	    computeNodeToNodeMappings(readNodes, readNodeCounts,
				      readPairs, cats, dubious, lengths,
				      false);
	computeNodeToNodeMappings(readNodes, readNodeCounts, readPairs,
				  cats, dubious, lengths, true);
	computeLocalNodeToNodeMappings();
	//computeLocalNodeToNodeMappings(readNodes, readNodeCounts, dubious);
	//fillUpLocalConnections(readNodes, readNodeCounts, lengths,
	//                     dubious);
	removeUnreliableConnections(cats);
	sortScaffold();

	// Clean up
	free(lengths);
	free(readNodeCounts);
	for (index = 0; index < sequenceCount(graph) + 1; index++)
		free(readNodes[index]);
	free(readNodes);
}

void cleanScaffoldMemory()
{
	free(scaffold);
	destroyRecycleBin(connectionMemory);
	connectionMemory = NULL;
}
