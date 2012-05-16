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

static PassageMarkerI nextMarker = NULL_IDX;
static PassageMarkerI nextTranscript = NULL_IDX;

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

////////////////////////////////////////
// Comparing PassageMarkers
////////////////////////////////////////

static int compareDirectedMarkers(PassageMarkerI A, PassageMarkerI B) {
	PassageMarkerI nextA = getNextInSequence(A);
	PassageMarkerI nextB = getNextInSequence(B);

	if (nextA && nextB) {
		if (getNode(nextA) != getNode(nextB))
			return 0;
		else
			return compareDirectedMarkers(nextA, nextB);
	} else if (nextA) 
		return 1;
	else if (nextB)
		return -1;
	else
		return 1;
}

static int compareMarkers(PassageMarkerI A, PassageMarkerI B) {
	int test1 = compareDirectedMarkers(A, B);
	int test2 = compareDirectedMarkers(getTwinMarker(A), getTwinMarker(B));

	if (test1 > 0 && test2 > 0)
		return 1;
	else if (test1 < 0 && test2 < 0)
		return -1;
	else
		return 0;
}

////////////////////////////////////////
// Removing PassageMarkers
////////////////////////////////////////

static void destroyTranscript(PassageMarkerI marker) {
	PassageMarkerI first, prev, ptr, next;

	first = marker;
	while((prev = getPreviousInSequence(first)))
		first = prev;

	for (ptr = first; ptr; ptr = next) {
		next = getNextInSequence(ptr);
		if (nextMarker == ptr || nextMarker == getTwinMarker(ptr))
			nextMarker = getNextInNode(nextMarker);
		if (nextTranscript == ptr || nextTranscript == getTwinMarker(ptr))
			nextTranscript = getNextInNode(nextTranscript);
		destroyPassageMarker(ptr);
	}
}

static boolean removeRedundantTranscripts_Marker(PassageMarkerI marker, boolean uniqueNodes) {
	PassageMarkerI marker2;
	if (getUniqueness(getNode(marker)) != uniqueNodes)
		return false;
	
	for (marker2 = getMarker(getNode(marker)); marker2; marker2 = nextMarker) {
		nextMarker = getNextInNode(marker2);

		if (uniqueNodes && getPassageMarkerStatus(marker2) == 1)
			continue;
		if (!uniqueNodes && getPassageMarkerStatus(marker2) == 2)
			continue;

		int comparison = compareMarkers(marker, marker2);

		if (comparison == 0)
			continue;
		else if (comparison > 0)
			destroyTranscript(marker2);
		else {
			destroyTranscript(marker);
			return true;
		}
	}	

	return false;
}

static void removeRedundantTranscripts_Transcript(PassageMarkerI marker, boolean uniqueNodes) {
	PassageMarkerI first, prev, ptr;

	if (getPassageMarkerStatus(marker))
		return; 

	first = marker;
	while((prev = getPreviousInSequence(first)))
		first = prev;

	if (uniqueNodes)
		for (ptr = first; ptr; ptr = getNextInSequence(ptr))
			setPassageMarkerStatus(ptr, 1);
	else
		for (ptr = first; ptr; ptr = getNextInSequence(ptr))
			setPassageMarkerStatus(ptr, 2);

	for (ptr = first; ptr; ptr = getNextInSequence(ptr))
		if (removeRedundantTranscripts_Marker(marker, uniqueNodes))
			break;
}

static void removeRedundantTranscripts_Node(Node * node, boolean uniqueNodes) {
	PassageMarkerI marker;
	
	if (getUniqueness(node) != uniqueNodes)
		return;

	for (marker = getMarker(node); marker; marker = nextTranscript) {
		nextTranscript = getNextInNode(marker);
		removeRedundantTranscripts_Transcript(marker, uniqueNodes);
	}
}

void removeRedundantTranscripts(Graph * argGraph) {
	IDnum index;
	graph = argGraph;
	
	velvetLog("Removing redundant transcripts\n");

	resetPassageMarkerStatuses(graph);
	defineUniqueNodes(graph);

	velvetLog("Scanning long nodes\n");
	for (index = 1; index <= nodeCount(graph); index++)
		removeRedundantTranscripts_Node(getNodeInGraph(graph, index), true);

	velvetLog("Scanning short nodes\n");
	for (index = 1; index <= nodeCount(graph); index++)
		removeRedundantTranscripts_Node(getNodeInGraph(graph, index), false);
}
