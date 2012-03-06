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
#ifndef _LOCUS_H_
#define _LOCUS_H_

#include "scaffold.h"

typedef struct locus_st Locus;
typedef struct transcript_st Transcript;
typedef struct event_st Event;

struct locus_st {
	IDnum contigCount;
	IDnum longContigCount;
	Node **contigs;
	Transcript *transcript;
	Event *event;
};

struct transcript_st {
	IDnum contigCount;
	Node **contigs;
	Coordinate *distances;
	double confidence;
	Transcript *next;
};

Locus *extractGraphLoci(Graph * graph, ReadSet * argReads,
			boolean * dubious, ShortLength * lengths,
			IDnum * locusCount, boolean scaffolding);

Connection *getReverseActiveConnection(Node * node);

// Modifying on the fly
void simplifyFromNode(Node * node, Locus * locus);
void removeNodeFromLocus(Node * node);
void renumberLocusNodes(Locus * locus);
void setLocusConnectionStatus(Locus * locus, boolean status);
void setLocusStatus(Locus * locus, boolean status);

// Constants
void setUnreliableConnectionCutoff_oases(int val);
void setDegreeCutoff(int val);
void setPairedThreshold(double pairedThreshold);

#endif
