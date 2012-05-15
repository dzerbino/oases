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
#ifndef _TRANSCRIPT_H_
#define _TRANSCRIPT_H_

typedef struct transcript_st Transcript;

#include "locus.h"

Transcript * newTranscript(IDnum contigCount, double confidence);
void produceTranscript(Locus * locus, IDnum nodesInList);
void addContigToTranscript(Transcript * transcript, Node * node, Coordinate distance);

// Getter/Setter
Transcript * getNextTranscript(Transcript * transcript);
void setNextTranscript(Transcript * transcript, Transcript * next);
IDnum getTranscriptContigCount(Transcript * transcript);
Node * getTranscriptContig(Transcript * transcript, IDnum index);
Coordinate getTranscriptDistance(Transcript * transcript, IDnum index);
double getConfidence(Transcript * transcript);

// Computation
void computeTranscripts(Graph * graph, Locus * loci, IDnum locusCount);

// Utility
void cleanTranscriptMemory(Locus * loci, IDnum locusCount);
void cleanLocusMemory(Locus * loci, IDnum locusCount);
#endif
