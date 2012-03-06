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
#ifndef _OASESEXPORT_H_
#define _OASESEXPORT_H_

void exportTranscripts(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength, Graph * graph);
void exportContigOrders(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength, Graph * graph);
void exportUnusedTranscriptReads(Graph* graph, Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char* directory);
IDnum usedTranscriptReads(Graph * graph, Coordinate minTransLength, Locus * loci, IDnum locusCount);
void exportAMOSTranscripts(Graph * graph, Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char * directory);
void exportTranscriptMappings(Locus * loci, IDnum locusCount, Graph * graph, ReadSet * reads, Coordinate minLength, char * directory);

#endif
