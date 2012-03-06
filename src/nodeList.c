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

static NodeList *markedNodes;
static RecycleBin *nodeListMemory = NULL;

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

NodeList *recordNode(Node * node)
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

boolean existsMarkedNode() {
	return markedNodes != NULL;
}

IDnum countMarkedNodes()
{
	NodeList *list;
	IDnum counter = 0;

	for (list = markedNodes; list != NULL; list = list->next)
		counter++;

	return counter;
}

Node *popNodeRecord()
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

void cleanNodeListMemory() {
	destroyRecycleBin(nodeListMemory);
	nodeListMemory = NULL;
}
