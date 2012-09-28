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

typedef struct locus_st Locus;

#include "scaffold.h"
#include "transcript.h"

Connection *getReverseActiveConnection(Node * node);

// Getters and setters
Transcript * getTranscript(Locus * locus);
void addTranscript(Locus * locus, Transcript * transcript);
IDnum getContigCount(Locus *locus);
void setContigCount(Locus *locus, IDnum count);
IDnum getLongContigCount(Locus *locus);
void setLongContigCount(Locus *locus, IDnum count);
Node * getContig(Locus * locus, IDnum index);
void addContig(Locus * locus, Node * contig);

// Modifying on the fly
void simplifyFromNode(Node * node, Locus * locus);
void removeNodeFromLocus(Node * node);
void renumberLocusNodes(Locus * locus);
void setLocusConnectionStatus(Locus * locus, boolean status);
void setLocusStatus(Locus * locus, boolean status);
void revert(Locus * locus);

// Constants
void setUnreliableConnectionCutoff_oases(int val);
void setDegreeCutoff(int val);
void setPairedThreshold(double pairedThreshold);

// Debugging
void exportLocusGraph(FILE * outfile, IDnum index, Locus * loci);

// Utility
Locus * allocateLocusArray(IDnum size);
Locus * reallocateLocusArray(Locus * array, IDnum newsize);
void clearLocus(Locus * locus);
void cleanLocusMemory(Locus * loci, IDnum locusCount);
Locus * getLocus(Locus * loci, IDnum index);
void printOasesConnections(Category * categories, char * filename);
#endif
