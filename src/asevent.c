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
	struct event_st *next;
};

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

