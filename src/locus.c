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
			clearLocus(locus);

			// Long contigs
			fillUpComponent(node);
			setLongContigCount(locus, countMarkedNodes());
			while (existsMarkedNode()) 
				addContig(locus, popNodeRecord());

			// Secondary contigs
			extendComponent(locus);
			setContigCount(locus, getLongContigCount(locus) + countMarkedNodes());
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

