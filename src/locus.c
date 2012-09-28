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
#include "concatenatedGraph.h"
#include "tightString.h"
#include "nodeList.h"
#include "locus.h"
#include "locallyCorrectedGraph2.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#define BLOCK_SIZE  100000
#define LENGTHCUTOFF 50

static Graph *graph = NULL;

struct locus_st {
	IDnum contigCount;
	IDnum longContigCount;
	Node **contigs;
	Transcript *transcript;
};

static void setNodeConnectionStatus(Node * node, boolean status)
{
	Connection *connect;

	for (connect = getConnection(node); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus
		    (getTwinNode(getConnectionDestination(connect))))
			setConnectionStatus(connect, status);
}

void setLocusConnectionStatus(Locus * locus, boolean status)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		setNodeConnectionStatus(locus->contigs[index], status);
}

void removeNodeFromLocus(Node * node) {
	setNodeStatus(node, false);
	setNodeConnectionStatus(node, false);
}

void renumberLocusNodes(Locus * locus) {
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

Connection *getReverseActiveConnection(Node * node)
{
	Connection *connect;

	for (connect = getConnection(getTwinNode(node)); connect;
	     connect = getNextConnection(connect))
		if (getNodeStatus(getConnectionDestination(connect)))
			return connect;

	return NULL;
}

void simplifyFromNode(Node * node, Locus * locus) {
	if (getReverseActiveConnection(node))
		return;

	correctGraphLocally2(node, locus);
}

void setLocusStatus(Locus * locus, boolean status)
{
	IDnum index;

	for (index = 0; index < locus->contigCount; index++)
		setSingleNodeStatus(locus->contigs[index], status);
}

void setDegreeCutoff(int val) {
	scaffold_setDegreeCutoff(val);
}

void setPairedThreshold(double pairedThreshold) {
	scaffold_setPairedThreshold(pairedThreshold);
}

Transcript * getTranscript(Locus * locus) {
	return locus->transcript;
}

void addTranscript(Locus * locus, Transcript * transcript) {
	setNextTranscript(transcript, locus->transcript);
	locus->transcript = transcript;
}

IDnum getContigCount(Locus *locus) {
	return locus->contigCount;
}

void setContigCount(Locus *locus, IDnum count) {
	if (locus->contigs)
		locus->contigs = reallocOrExit(locus->contigs, count, Node *);
	else
		locus->contigs = callocOrExit(count, Node *);
}

IDnum getLongContigCount(Locus *locus) {
	return locus->longContigCount;
}

void setLongContigCount(Locus *locus, IDnum count) {
	locus->longContigCount = count;
	locus->contigs = callocOrExit(count, Node *);
}

Node * getContig(Locus * locus, IDnum index) {
	return locus->contigs[index];
}

void addContig(Locus * locus, Node * contig) {
	locus->contigs[locus->contigCount++] = contig;
}

void cleanLocusMemory(Locus * loci, IDnum locusCount)
{
	IDnum index;

	cleanTranscriptMemory(loci, locusCount);
	for (index = 0; index < locusCount; index++)
		free(loci[index].contigs);
	free(loci);
	cleanScaffoldMemory();
	cleanNodeListMemory();
}

void clearLocus(Locus * locus) {
	locus->contigs = NULL;
	locus->contigCount = 0;
	locus->longContigCount = 0;
	locus->transcript = NULL;
}

Locus * getLocus(Locus * loci, IDnum index) {
	return loci + index;
}

Locus * allocateLocusArray(IDnum size) {
	return callocOrExit(size, Locus);
}

Locus * reallocateLocusArray(Locus * locus, IDnum size) {
	return reallocOrExit(locus, size, Locus);
}

void exportLocusGraph(FILE * file, IDnum index, Locus * loci)
{
	Locus *locus = getLocus(loci, index);
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
	for (i = 0; i < getLongContigCount(locus); i++) {
		nodeID = (long long) getNodeID(getContig(locus, i));
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
	for (i = getLongContigCount(locus); i < getContigCount(locus); i++) {
		nodeID = (long long) getNodeID(getContig(locus, i));
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
	for (i = 0; i < getContigCount(locus); i++) {
		node = getContig(locus, i);
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

void revert(Locus * locus)
{
	IDnum index;

	for (index = 0; index < getContigCount(locus); index++)
		locus->contigs[index] = getTwinNode(getContig(locus, index));
}

void printOasesConnections(Category * categories, char * filename) {
	FILE * file = fopen(filename, "w");
	if (!file) 
		exitErrorf(-1, false, "Could not open %s!i Exiting...\n", filename);
	printOasesConnections2(categories, file);
}
