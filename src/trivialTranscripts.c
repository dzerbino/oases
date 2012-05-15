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

/////////////////////////////////////////
// Locus degree distribution
/////////////////////////////////////////
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

void getDegreeDistribution(Locus * locus, IDnum * distribution)
{
	IDnum index;

	for (index = 0; index < 4; index++)
		distribution[index] = 0;

	for (index = 0; index < getContigCount(locus); index++)
		updateDistribution(distribution, getContig(locus, index));
}

#define UNKNOWN 0
#define LINEAR 1
#define FORK 2
#define BUBBLE 3

boolean hasPlausibleTranscripts(IDnum * d)
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

///////////////////////////////////////////
// Extracting trivial transcripts
///////////////////////////////////////////
IDnum extractLinearPath(Node * maxNode)
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

	for (index = 0; index < getContigCount(locus); index++)
		if (hasNoActiveConnections(getContig(locus, index)))
			return getContig(locus, index);

	return NULL;
}

void addLinearTranscriptToLocus(Locus * locus)
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

void addBubbleTranscriptsToLocus(Locus * locus)
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

	for (index = 0; index < getContigCount(locus); index++) {
		if (hasNoActiveConnections(getContig(locus, index))) {
			if (firstFound)
				return getContig(locus, index);
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

void addForkedTranscriptsToLocus(Locus * locus)
{
	if (onlyOneEndPoint(locus))
		addBubbleTranscriptsToLocus(locus);
	else
		add3primeForkedTranscriptsToLocus(locus);
}
