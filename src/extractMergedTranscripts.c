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

#include "globals.h"
#include "utility.h"
#include "graph.h"
#include "passageMarker.h"
#include "transcript.h"
#include "recycleBin.h"
#include "locallyCorrectedGraph.h"
#include "transcript.h"
#include "scaffold.h"
#include "nodeList.h"

static Graph * graph = NULL;

////////////////////////////////////////
// PassageMarker status 
////////////////////////////////////////

static void resetPassageMarkerStatuses_Node(Node * node) {
	PassageMarkerI marker;

	for (marker = getMarker(node); marker; marker = getNextInNode(marker))
		setPassageMarkerStatus(marker, false);
}

static void resetPassageMarkerStatuses(Graph * graph) {
	IDnum index;

	for (index = 1; index <= nodeCount(graph); index++)
		resetPassageMarkerStatuses_Node(getNodeInGraph(graph, index));
}

/////////////////////////////////////////
// Counting loci
////////////////////////////////////////

static Node * getNextDestination(PassageMarkerI marker) {
	PassageMarkerI ptr;

	for (ptr = getNextInSequence(marker); ptr; ptr = getNextInSequence(ptr))
		if (getUniqueness(getNode(ptr)))
			return getNode(ptr);
	return NULL;
}

static void propagateComponent(Node * node)
{
	PassageMarkerI marker;

	if (!node || getNodeStatus(node) || !getUniqueness(node))
		return;

	setNodeStatus(node, true);

	for (marker = getMarker(node); marker; marker = getNextInNode(marker))
		propagateComponent(getNextDestination(marker));

	for (marker = getMarker(getTwinNode(node)); marker; marker = getNextInNode(marker))
		propagateComponent(getNextDestination(marker));
}

static IDnum countConnectedComponents(Graph * argGraph)
{
	IDnum index;
	IDnum count = 0;
	Node *node;

	graph = argGraph;

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

/////////////////////////////////////////
// Creating loci
////////////////////////////////////////

static void extendComponentToNode(Node * node)
{
	if (!node || getNodeStatus(node))
		return;

	setSingleNodeStatus(node, true);
	recordNode(node);
}

static void extendComponentFromNode(Node * node)
{
	PassageMarkerI marker;
	Node * next;

	for (marker = getMarker(node); marker; marker = getNextInNode(marker))
		extendComponentToNode(getNextDestination(marker));

	for (marker = getMarker(getTwinNode(node)); marker; marker = getNextInNode(marker))
		if ((next = getNextDestination(marker)))
			extendComponentToNode(getTwinNode(next));
}

static void extendComponent(Locus * locus)
{
	IDnum index;

	for (index = 0; index < getLongContigCount(locus); index++)
		extendComponentFromNode(getContig(locus, index));
}

static void fillUpComponent(Node * node)
{
	PassageMarkerI marker;
	Node * next;

	if (!node || getNodeStatus(node) || !getUniqueness(node))
		return;
	setSingleNodeStatus(node, true);
	recordNode(node);

	for (marker = getMarker(node); marker; marker = getNextInNode(marker))
		fillUpComponent(getNextDestination(marker));

	for (marker = getMarker(getTwinNode(node)); marker; marker = getNextInNode(marker))
		if ((next = getNextDestination(marker)))
			fillUpComponent(getTwinNode(next));
}

static Locus *extractConnectedComponents(IDnum locusCount)
{
	Locus *loci = allocateLocusArray(locusCount);
	Locus *locus;
	IDnum index;
	IDnum locusIndex = 0;
	IDnum nodeIndex;
	Node *node;

	resetNodeStatus(graph);

	for (index = 1; index <= nodeCount(graph); index++) {
		node = getNodeInGraph(graph, index);
		if (!getNodeStatus(node) && getUniqueness(node)) {
			locus = getLocus(loci, locusIndex++);

			// Long contigs
			fillUpComponent(node);
			setLongContigCount(locus, countMarkedNodes());
			while (existsMarkedNode()) 
				addContig(locus, popNodeRecord());

			// Secondary contigs
			extendComponent(locus);
			setContigCount(locus, getLongContigCount(locus) + countMarkedNodes());
			// DEBUG
			if (getLongContigCount(locus) > 1 && countMarkedNodes() == 0)
				abort();
			while (existsMarkedNode())
				addContig(locus, popNodeRecord());

			// Mark primary nodes so that their twins are not reused
			for (nodeIndex = 0;
			     nodeIndex < getLongContigCount(locus);
			     nodeIndex++)
				setNodeStatus(getContig(locus, nodeIndex), true);

			// Unmark secondary nodes so that they are available to other loci
			for (nodeIndex = getLongContigCount(locus);
			     nodeIndex < getContigCount(locus); nodeIndex++)
				setNodeStatus(getContig(locus, nodeIndex), false);
		}
	}

	return loci;
}

Locus *reextractGraphLoci(Graph * graph, IDnum * locusCount) {
	velvetLog("Re-extracting loci from connection graph...\n");

	*locusCount = countConnectedComponents(graph);
	velvetLog("Counted %li mRNA loci\n", (long) *locusCount);

	return extractConnectedComponents(*locusCount);
}

/////////////////////////////////////////
// Extract transcripts
/////////////////////////////////////////

static IDnum listTranscriptNodes(PassageMarkerI marker) {
	PassageMarkerI current;
	IDnum nodesInList = 0;

	for (current = marker; current; current = getPreviousInSequence(current)) {
		recordNode(getNode(current));
		nodesInList++;
	}

	return nodesInList;
}

static void produceTranscript_Merge(Locus * locus, IDnum nodesInList)
{
	Node *node;

	Transcript *transcript = newTranscript(nodesInList, ((double) nodesInList) / getContigCount(locus));

	while ((node = popNodeRecord())) 
		addContigToTranscript(transcript, node, 0);

	addTranscript(locus, transcript);
}

static void addTranscriptToLocus_Marker(Locus * locus, PassageMarkerI marker) {
	PassageMarkerI last, next, ptr;

	if (getPassageMarkerStatus(marker))
		return;
	setPassageMarkerStatus(marker, true);

	last = marker;
	while((next = getNextInSequence(last)))
		last = next;

	for (ptr = last; ptr; ptr = getPreviousInSequence(ptr))
		setPassageMarkerStatus(ptr, true);

	produceTranscript_Merge(locus, listTranscriptNodes(last));
}

static void addTranscriptToLocus_Node(Locus * locus, IDnum index) {
	PassageMarkerI marker;

	for (marker = getMarker(getContig(locus, index)); marker; marker = getNextInNode(marker))
		addTranscriptToLocus_Marker(locus, marker);
}

static void addTranscriptToLocus(Locus * locus)
{
	IDnum index;

	for (index = 0; index < getContigCount(locus); index++)
		addTranscriptToLocus_Node(locus, index);
}

static int compareNodes(const void * A, const void * B) {
	Node * a = *((Node **) A);
	Node * b = *((Node **) B);
	if (getNodeID(a) > getNodeID(b))
		return 1;
	else if (getNodeID(a) < getNodeID(b))
		return -1;
	else
		return 0;
}

static void addOrphanedTranscriptToLoci(Locus * locus, PassageMarkerI marker) {
	IDnum index;
	PassageMarkerI current, start;
	Node ** nodes;
	IDnum length = 0;
	IDnum counter = 0;

	clearLocus(locus);

	// Go to start of sequence
	for (start = marker; getPreviousInSequence(start); start = getPreviousInSequence(start))
		;

	// Count nodes
	for (current = start; current; current = getNextInSequence(current))
		counter++;

	// Copy nodes
	nodes = callocOrExit(counter, Node *);
	for (current = start; current; current = getNextInSequence(current))
		nodes[length++] = getNode(current);	

	// Sort nodes
	qsort(nodes, length, sizeof(Node *), compareNodes);
	
	// Count single occurences
	counter = 1;
	for (index = 1; index < length; index++)
		if (nodes[index] != nodes[index - 1])
			counter++;	

	// Fill nodes array
	setContigCount(locus, counter);
	addContig(locus, nodes[0]);
	for (index = 1; index < length; index++)
		if (nodes[index] != nodes[index - 1])
			addContig(locus, nodes[index]);

	addTranscriptToLocus_Marker(locus, marker);
	free(nodes);
}

void recomputeTranscripts(Locus ** loci, IDnum * locusCount) {
	IDnum index;
	Locus *locus;
	IDnum counter = 0;
	PassageMarkerI marker;
	boolean * markedTranscripts;

	velvetLog("Computing merged transcripts\n");

	resetPassageMarkerStatuses(graph);

	for (index = 0; index < *locusCount; index++) {
		locus = getLocus(*loci, index);
		addTranscriptToLocus(locus);
	}

	velvetLog("Catching orphaned transcripts\n");

	markedTranscripts = callocOrExit(sequenceCount(graph) + 1, boolean);
	
	for (index = 1; index <= nodeCount(graph); index++)
		for (marker = getMarker(getNodeInGraph(graph, index)); marker; marker = getNextInNode(marker))
			if (!getPassageMarkerStatus(marker))
				markedTranscripts[getAbsolutePassMarkerSeqID(marker)] = true;

	for (index = 1; index <= sequenceCount(graph); index++) {
		if (markedTranscripts[index]) {
			markedTranscripts[index] = false;
			counter++;
		}
	}

	velvetLog("Found %li missing transcripts\n", (long) counter);

	*loci = reallocateLocusArray(*loci, *locusCount + counter);

	for (index = 1; index <= nodeCount(graph); index++) {
		for (marker = getMarker(getNodeInGraph(graph, index)); marker; marker = getNextInNode(marker)) {
			if (!getPassageMarkerStatus(marker) && !markedTranscripts[getAbsolutePassMarkerSeqID(marker)]) {
				addOrphanedTranscriptToLoci( getLocus(*loci, (*locusCount)++), marker);
				markedTranscripts[getAbsolutePassMarkerSeqID(marker)] = true;
			}
		}
	}

	free(markedTranscripts);
}
