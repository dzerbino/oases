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
#include "nodeList.h"
#include "locus.h"
#include "transcript.h"

#define ADENINE 0
#define CYTOSINE 1
#define GUANINE 2
#define THYMINE 3

#define BLOCK_SIZE  100000
#define LENGTHCUTOFF 50

static Graph * graph = NULL;

static IDnum abs_id(IDnum id)
{
	return id > 0 ? id : -id;
}

static void exportNodeSequence(Node * node, FILE * outfile)
{
	int wordShift = getWordLength(graph) - 1;
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	Coordinate length = getNodeLength(node) + wordShift;
	Coordinate start;
	char str[100];

	for (start = 0; start <= length; start += 60) {
		strncpy(str, &(string[start]), 60);
		str[60] = '\0';
		fprintf(outfile, "%s\n", str);
	}
	free(string);
}

static void exportLocusNode(IDnum index, Node * node, FILE * outfile)
{
	fprintf(outfile, ">Locus_%li_Node_%li\n", (long) index + 1,
		(long) abs_id(getNodeID(node)));
	exportNodeSequence(node, outfile);
}

static void exportLocusNodes(IDnum locusIndex, Locus * locus,
			     FILE * outfile)
{
	IDnum index;

	for (index = 0; index < getContigCount(locus); index++)
		exportLocusNode(locusIndex, getContig(locus, index),
				outfile);
}

static Coordinate getTranscriptLength(Transcript * transcript) {
	IDnum index;
	Coordinate totalLength = 0;

	if (getTranscriptContigCount(transcript) == 0)
		return 0;

	totalLength = getWordLength(graph) - 1;

	for (index = 0; index < getTranscriptContigCount(transcript); index++) {
		totalLength += getNodeLength(getTranscriptContig(transcript, index));
		if (index < getTranscriptContigCount(transcript) - 1) {
		        // Ignoring NCBI request for 10bp gaps
			//if (getTranscriptDistance(transcript, index) < getWordLength(graph) && getNodeLength(getTranscriptContig(transcript, index+1)) >= getWordLength(graph) - 1) {
			//	totalLength += getTranscriptDistance(transcript, index);
			//} else if (getTranscriptDistance(transcript, index) > 0 && getTranscriptDistance(transcript, index) < 10) {
			//	totalLength += 10;
			//} else {
			totalLength += getTranscriptDistance(transcript, index);
			//}
		}
	}

	return totalLength;
}


static void exportTranscriptContigs(Transcript * transcript, IDnum locusID,
				    IDnum transID, IDnum transcriptCount, FILE * outfile)
{
	IDnum index;
	Coordinate totalLength = 0;

	if (getTranscriptContigCount(transcript) == 0)
		return;

	// Header
	fprintf(outfile, ">Locus_%li_Transcript_%li/%li_Confidence_%.3f_Length_%li\n",
		(long) locusID + 1, (long) transID + 1, (long) transcriptCount, getConfidence(transcript), (long) getTranscriptLength(transcript));

	totalLength = getWordLength(graph) - 1;

	// Sequence
	for (index = 0; index < getTranscriptContigCount(transcript); index++) {
		totalLength += getNodeLength(getTranscriptContig(transcript, index));
		fprintf(outfile, "%li:%lli",
			(long) getNodeID(getTranscriptContig(transcript, index)),
			(long long) totalLength);
		if (index < getTranscriptContigCount(transcript) - 1) {
			if (getTranscriptDistance(transcript, index) < getWordLength(graph) && getNodeLength(getTranscriptContig(transcript, index+1)) >= getWordLength(graph) - 1) {
				fprintf(outfile, "-(0)->");
				totalLength += getTranscriptDistance(transcript, index);
			// Ignoring NCBI request for 10bp gaps
			//} else if (getTranscriptDistance(transcript, index) > 0 && getTranscriptDistance(transcript, index) < 10) {
			//	fprintf(outfile, "-(10)->");
			//	totalLength += 10;
			} else {
				fprintf(outfile, "-(%li)->",
					(long) getTranscriptDistance(transcript, index));
				totalLength += getTranscriptDistance(transcript, index);
			}
		}
	}

	fprintf(outfile, "\n");
}

static void exportLocusContigs(IDnum locusID, Locus * locus,
			       FILE * outfile, Coordinate minTransLength)
{
	IDnum index = 0;
	Transcript *transcript;
	IDnum transcriptCount = 0;

	for (transcript = getTranscript(locus); transcript != NULL;
	     transcript = getNextTranscript(transcript))
		if (getTranscriptLength(transcript) >= minTransLength)
			transcriptCount++;

	exportLocusNodes(locusID, locus, outfile);
	for (transcript = getTranscript(locus); transcript != NULL;
	     transcript = getNextTranscript(transcript))
		if (getTranscriptLength(transcript) >= minTransLength)
			exportTranscriptContigs(transcript, locusID, index++, transcriptCount,
						outfile);
}

void exportContigOrders(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength, Graph * argGraph)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;
	graph = argGraph;

	if (outfile)
		velvetLog("Exporting transcript contigs to %s\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	for (index = 0; index < locusCount; index++)
		exportLocusContigs(index, getLocus(loci, index), outfile, minTransLength);

	fclose(outfile);

}

ReadSet *importEmptyReadSet(char *filename, Coordinate ** lengthsPtr,
			    int wordLength)
{
	FILE *file = fopen(filename, "r");
	const int maxline = 5000;
	char line[5000];
	IDnum sequenceCount, sequenceIndex;
	ReadSet *reads;
	short int temp_short;
	int lengthOffset = wordLength - 1;
	Coordinate bpCount = 0;

	if (file != NULL)
		velvetLog("Reading read set file %s;\n", filename);
	else
		exitErrorf(EXIT_FAILURE, true, "Could not open %s",
			   filename);

	reads = newReadSet();

	// Count number of separate sequences
	sequenceCount = 0;
	while (fgets(line, maxline, file) != NULL)
		if (line[0] == '>')
			sequenceCount++;
	fclose(file);
	velvetLog("%ld sequences found\n", (long) sequenceCount);

	reads->readCount = sequenceCount;

	if (reads->readCount == 0) {
		reads->sequences = NULL;
		reads->categories = NULL;
		return reads;
	}

	*lengthsPtr = callocOrExit(sequenceCount, Coordinate);
	reads->categories = callocOrExit(sequenceCount, Category);
	// Counting base pair length of each sequence:
	file = fopen(filename, "r");
	sequenceIndex = -1;
	while (fgets(line, maxline, file) != NULL) {
		if (line[0] == '>') {

			// Reading category info
			sscanf(line, "%*[^\t]\t%*[^\t]\t%hd", &temp_short);
			reads->categories[sequenceIndex + 1] =
			    (Category) temp_short;
			if (sequenceIndex != -1)
				(*lengthsPtr)[sequenceIndex] =
				    bpCount - lengthOffset;
			sequenceIndex++;
			bpCount = 0;
		} else {
			bpCount += (Coordinate) strlen(line) - 1;
		}
	}
	(*lengthsPtr)[sequenceIndex] = bpCount - lengthOffset;
	fclose(file);

	puts("Done");
	return reads;

}

static void exportAMOSLib(FILE * outfile, Graph * graph, Category cat)
{
	Coordinate distance = getInsertLength(graph, cat * 2);
	double variance = getInsertLength_var(graph, cat * 2);

	if (distance == -1)
		return;

	fprintf(outfile, "{LIB\n");
	fprintf(outfile, "iid:%d\n", (int) (cat + 1));
	fprintf(outfile, "{DST\n");
	fprintf(outfile, "mea:%lld\n", (long long) distance);
	fprintf(outfile, "std:%lld\n", (long long) sqrt(variance));
	fprintf(outfile, "}\n");
	fprintf(outfile, "}\n");
}

static void exportAMOSMarker(FILE * outfile, PassageMarkerI marker,
			     Coordinate nodeLength, Coordinate offset,
			     int wordShift, boolean firstUniqueMet)
{
	Coordinate sequenceStart, sequenceFinish;

	sequenceStart = getPassageMarkerStart(marker);
	sequenceFinish = getPassageMarkerFinish(marker);

	if (firstUniqueMet == -1) {
		if (getPassageMarkerSequenceID(marker) < 0) {
			sequenceFinish += wordShift;
			sequenceStart += wordShift;
		}
	} else if (firstUniqueMet == 0) {
		if (getPassageMarkerSequenceID(marker) > 0) 
			sequenceFinish += wordShift;
		else 
			sequenceStart += wordShift;
	} else if (firstUniqueMet == 1) {
		if (getPassageMarkerSequenceID(marker) > 0) {
			sequenceFinish += wordShift;
			sequenceStart += wordShift;
		}
	}

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%ld\n", (long) getAbsolutePassMarkerSeqID(marker));
	fprintf(outfile, "off:%lld\n", (long long) (offset + getStartOffset(marker)));
	fprintf(outfile, "clr:%lld,%lld\n", (long long) sequenceStart, (long long) sequenceFinish);
	fprintf(outfile, "}\n");
}

static void exportAMOSShortMarker(FILE * outfile, ShortReadMarker * marker,
				  ReadSet * reads, Coordinate offset, boolean firstUniqueMet, int wordShift)
{
	Coordinate read_offset =
	    getShortReadMarkerPosition(marker) -
	    getShortReadMarkerOffset(marker)
	    + offset;
	TightString *sequence = getTightStringInArray(reads->tSequences, getShortReadMarkerID(marker) - 1);

	if (firstUniqueMet == 1)
		read_offset -= wordShift;

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%ld\n", (long) getShortReadMarkerID(marker));
	fprintf(outfile, "off:%lld\n", (long long) read_offset);
	fprintf(outfile, "clr:0,%lld\n", (long long) getLength(sequence));
	fprintf(outfile, "}\n");
}

static void exportAMOSReverseShortMarker(FILE * outfile,
					 ShortReadMarker * marker,
					 Coordinate nodeLength,
					 int wordShift, ReadSet * reads,
					 Coordinate offset, boolean firstUniqueMet)
{
	TightString *sequence = getTightStringInArray(reads->tSequences, getShortReadMarkerID(marker) - 1);

	Coordinate read_offset =
	    nodeLength - getShortReadMarkerPosition(marker) +
	    getShortReadMarkerOffset(marker) - getLength(sequence) +
	    wordShift + offset;

	if (firstUniqueMet == 1)
		read_offset -= wordShift;

	if (getShortReadMarkerPosition(marker) == -1)
		return;

	fprintf(outfile, "{TLE\n");
	fprintf(outfile, "src:%ld\n", (long) getShortReadMarkerID(marker));
	fprintf(outfile, "off:%lld\n", (long long) read_offset);
	fprintf(outfile, "clr:%lld,0\n", (long long) getLength(sequence));
	fprintf(outfile, "}\n");
}

static char * nSequence(Coordinate length) {
	char * sequence = callocOrExit(length + 1, char);
	Coordinate position;

	for (position = 0; position < length; position++) 
		sequence[position] = 'N';
	sequence[length] = '\0';

	return sequence;
}

static char revComp(char base) {
	if (base == 'A')
		return 'T';
	if (base == 'T')
		return 'A';
	if (base == 'G')
		return 'C';
	if (base == 'C')
		return 'G';
	return 'N';
}

static void printFastA(FILE * outfile, char * sequence) {
	Coordinate position;
	size_t length = strlen(sequence);
	char str[100];

	for (position = 0; position < length; position += strlen(str)) {
		strncpy(str, sequence + position, 60);
		str[60] = '\0';
		fprintf(outfile, "%s\n", str);
	}
}

static void addLongNodeToSequence(char * sequence, Node * node, Coordinate offset) {
	Coordinate position;
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	int wordShift = getWordLength(graph) - 1;

	for (position = 0; position < strlen(string); position++) 
	    if (sequence[offset + position - wordShift] == 'N')
		sequence[offset + position - wordShift] = string[position]; 
	free(string);
}

static void addShortNodeToSequence(char * sequence, Node * node, Coordinate offset) {
	Coordinate position;
	char *string = expandNodeFragment(node, 0, getNodeLength(node),
					  getWordLength(graph));
	int wordShift = getWordLength(graph) - 1;

	for (position = 0; position < strlen(string); position++) 
	    if (sequence[offset + position] == 'N')
		sequence[offset + position] = string[position]; 

	free(string);
	string = expandNodeFragment(getTwinNode(node), 0, getNodeLength(node),
				      getWordLength(graph));

	for (position = 0; position < strlen(string); position++) 
	    if (sequence[offset - wordShift + getNodeLength(node) - 1 - position] == 'N')
		sequence[offset - wordShift + getNodeLength(node) - 1 - position] = revComp(string[position]); 
	free(string);
}

static void addNodeToSequence(char * sequence, Node * node, Coordinate offset) {
	int wordShift = getWordLength(graph) - 1;

	// Add sequence
	if (getNodeLength(node) >= wordShift) 
		addLongNodeToSequence(sequence, node, offset);
	else
		addShortNodeToSequence(sequence, node, offset);
}

static void printContigSequence(FILE * outfile, Transcript * transcript, IDnum startIndex, IDnum finishIndex, Coordinate * offsets, Coordinate * contigLengths) 
{
	IDnum index;
	Coordinate offset = getWordLength(graph) - 1;
	char * sequence = nSequence(getTranscriptLength(transcript)); 

	// Extract sequence
	for (index = startIndex; index <= finishIndex; index++) {
	    addNodeToSequence(sequence, getTranscriptContig(transcript, index), offset);
	    offset += getNodeLength(getTranscriptContig(transcript, index));
	    if (index < getTranscriptContigCount(transcript) - 1)
		    offset += getTranscriptDistance(transcript, index);
	}
	sequence[offset] = '\0';

	// Count initial N's
	for (offset = 0; offset < strlen(sequence); offset++)
		if (sequence[offset] != 'N')
			break;
	offsets[finishIndex] += offset;
	contigLengths[finishIndex] = strlen(sequence) - offset;
	
	printFastA(outfile, sequence + offset);
	free(sequence);
}

static void exportAMOSContig(FILE * outfile, ReadSet * reads, Transcript * transcript, 
			     IDnum startIndex, IDnum finishIndex, Graph * graph, 
			     IDnum locusID, IDnum transcriptID, Coordinate * offsets, Coordinate * contigLengths, IDnum internalIndex, IDnum iid)
{
	Coordinate start;
	PassageMarkerI marker;
	ShortReadMarker *shortMarkerArray, *shortMarker;
	Coordinate index, maxIndex;
	int wordShift = getWordLength(graph) - 1;
	Coordinate offset = 0;
	IDnum nodeIndex;
	Node * node;
	boolean firstUniqueMet = false;
	int column = 0;
	Coordinate totalLength = 0;

	// Check if any material to work with
	for (nodeIndex = startIndex; nodeIndex <= finishIndex; nodeIndex++)
		totalLength += getNodeLength(getTranscriptContig(transcript, nodeIndex));
	if (totalLength == 0)
		return;

	// Compute offset within transcript
	for (nodeIndex = 0; nodeIndex < startIndex; nodeIndex++)
		offsets[finishIndex] += getNodeLength(getTranscriptContig(transcript, nodeIndex)) + getTranscriptDistance(transcript, nodeIndex);

	fprintf(outfile, "{CTG\n");
	fprintf(outfile, "iid:%ld\n", (long) iid);
	fprintf(outfile, "eid:%ld-%ld-%ld\n", (long) locusID, (long) transcriptID, (long) internalIndex);

	fprintf(outfile, "seq:\n");
	printContigSequence(outfile, transcript, startIndex, finishIndex, offsets, contigLengths);

	fprintf(outfile, "qlt:\n");
	column = 0;
	for (start = 0; start < contigLengths[finishIndex]; start++) {
		fprintf(outfile, "I");

		if (column++ == 60) {
			fprintf(outfile, "\n");
			column = 0;
		}
	}
	fprintf(outfile, "\n.\n");

	firstUniqueMet = -1;
	for (nodeIndex = startIndex; nodeIndex <= finishIndex; nodeIndex++) {
		node = getTranscriptContig(transcript, nodeIndex);
		if (firstUniqueMet == 0)
			firstUniqueMet = 1;
		else if (firstUniqueMet == -1 && getNodeLength(node) >= wordShift) 
			firstUniqueMet = 0;

		for (marker = getMarker(node); marker != NULL_IDX;
		     marker = getNextInNode(marker))
			exportAMOSMarker(outfile, marker, getNodeLength(node),
					 offset, wordShift, firstUniqueMet);
		

		if (readStartsAreActivated(graph)) {
			shortMarkerArray = getNodeReads(node, graph);
			maxIndex = getNodeReadCount(node, graph);
			for (index = 0; index < maxIndex; index++) {
				shortMarker =
				    getShortReadMarkerAtIndex(shortMarkerArray,
							      index);
				exportAMOSShortMarker(outfile, shortMarker, reads,
						      offset, firstUniqueMet, wordShift);
			}

			shortMarkerArray = getNodeReads(getTwinNode(node), graph);
			maxIndex = getNodeReadCount(getTwinNode(node), graph);
			for (index = 0; index < maxIndex; index++) {
				shortMarker =
				    getShortReadMarkerAtIndex(shortMarkerArray,
							      index);
				exportAMOSReverseShortMarker(outfile, shortMarker,
							     getNodeLength(node),
							     wordShift, reads,
							     offset, firstUniqueMet);
			}
		}

		offset += getNodeLength(node);
	}

	fprintf(outfile, "}\n");
}

static void exportAMOSTranscript(FILE* outfile, ReadSet * reads, Transcript * transcript, IDnum locusID, IDnum transcriptID, Graph * graph) 
{
	IDnum smallIndex = 0;
	static IDnum iid = 1;
	IDnum contigIndex = iid;
	IDnum nodeIndex;
	IDnum startIndex = 0;
	Coordinate * contigLengths = callocOrExit(getTranscriptContigCount(transcript), Coordinate);
	Coordinate * offsets = callocOrExit(getTranscriptContigCount(transcript), Coordinate);

	for (nodeIndex = 0; nodeIndex < getTranscriptContigCount(transcript); nodeIndex++)
		if (nodeIndex == getTranscriptContigCount(transcript) - 1 || getTranscriptDistance(transcript, nodeIndex) > 0)
			    exportAMOSContig(outfile, reads, transcript, startIndex, nodeIndex, graph, 
						     locusID, transcriptID, offsets, contigLengths, smallIndex++, iid++);

	fprintf(outfile, "{SCF\n");
	fprintf(outfile, "eid:%ld-%ld\n", (long) locusID, (long) transcriptID);

	for (nodeIndex = 0; nodeIndex < getTranscriptContigCount(transcript); nodeIndex++) {
		if (nodeIndex == getTranscriptContigCount(transcript) - 1 || getTranscriptDistance(transcript, nodeIndex) > 0) {
			if (contigLengths[nodeIndex] > 0) {
				fprintf(outfile, "{TLE\n");
				fprintf(outfile, "off:%lld\n", (long long) offsets[nodeIndex]);
				fprintf(outfile, "clr:0,%lld\n", (long long) contigLengths[nodeIndex]);
				fprintf(outfile, "src:%ld\n", (long) contigIndex++);
				fprintf(outfile, "}\n");
			}
		} 
	}

	fprintf(outfile, "}\n");
	free(offsets);
	free(contigLengths);
}

static void exportAMOSLocus(FILE * outfile, ReadSet * reads, Locus * locus, IDnum locusID, Coordinate minTransLength,
			   Graph * graph)
{
	Transcript * transcript;
	IDnum transcriptID = 0;

	for (transcript = getTranscript(locus); transcript; transcript = getNextTranscript(transcript)) 
		if (getTranscriptLength(transcript) >= minTransLength)
			exportAMOSTranscript(outfile, reads, transcript, locusID, transcriptID++, graph);
}

static void exportAMOSRead(FILE * outfile, TightString * tString,
			   IDnum index, IDnum frg_index)
{
	Coordinate start, finish;
	char str[100];

	fprintf(outfile, "{RED\n");
	fprintf(outfile, "iid:%ld\n", (long) index);
	fprintf(outfile, "eid:%ld\n", (long) index);
	if (frg_index > 0)
		fprintf(outfile, "frg:%ld\n", (long) frg_index);

	fprintf(outfile, "seq:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "qlt:\n");
	start = 0;
	while (start <= getLength(tString)) {
		finish = start + 60;
		readTightStringFragment(tString, start, finish, str);
		fprintf(outfile, "%s\n", str);
		start = finish;
	}
	fprintf(outfile, ".\n");

	fprintf(outfile, "}\n");
}

void exportAMOSTranscripts(Graph * graph,
		       Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char * directory)
{
	IDnum index;
	Category cat;
	Locus * locus;
	FILE *outfile;
	char filename[10000];

	strcpy(filename, directory);
	strcat(filename, "/oases_asm.afg");
	velvetLog("Writing into AMOS file %s...\n", filename);
	outfile = fopen(filename, "w");

	if (outfile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to AMOS file %s",
		       filename);

	for (cat = 0; cat <= CATEGORIES; cat++)
		exportAMOSLib(outfile, graph, cat);

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			fprintf(outfile, "{FRG\n");
			fprintf(outfile, "lib:%d\n",
				(int) ((reads->categories[index - 1] / 2) + 1));
			fprintf(outfile, "rds:%ld,%ld\n", (long) index,
				(long) index + 1);
			fprintf(outfile, "eid:%ld\n", (long) index);
			fprintf(outfile, "iid:%ld\n", (long) index);
			fprintf(outfile, "typ:I\n");
			fprintf(outfile, "}\n");
			index++;
		}
	}

	for (index = 1; index <= reads->readCount; index++) {
		if (reads->categories[index - 1] % 2 != 0 &&
		    getInsertLength(graph,
				    reads->categories[index - 1]) >= 0) {
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       index);
			index++;
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       index - 1);
		} else {
			exportAMOSRead(outfile,
				       getTightStringInArray(reads->tSequences, index - 1), index,
				       -1);
		}
	}

	for (index = 0; index < locusCount; index++) {
		locus = getLocus(loci, index);

		if (locus == NULL)
			continue;

		exportAMOSLocus(outfile, reads, locus, index, minTransLength, graph);
	}

	fclose(outfile);

}

static void markUsedReads(Node * node, boolean * used)
{
	IDnum readID;
	ShortReadMarker * shortReadArray, * shortReadMarker;
	IDnum shortReadCount, shortReadIndex;
	PassageMarkerI marker;

	if (node == NULL || getNodeStatus(node))
		return;
	else
		setNodeStatus(node, true);
	
	// Long reads
	for(marker = getMarker(node); marker != NULL_IDX; marker = getNextInNode(marker)) {
		readID = getPassageMarkerSequenceID(marker);
		if (readID < 0)
			readID = -readID;
		used[readID] = true;	
	}	

	// Short reads		
	if (!readStartsAreActivated(graph))
		return;

	shortReadArray = getNodeReads(node, graph);
	shortReadCount = getNodeReadCount(node, graph);
	for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
		shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
		readID = getShortReadMarkerID(shortReadMarker);
		used[readID] = true;	
	}
	
	shortReadArray = getNodeReads(getTwinNode(node), graph);
	shortReadCount = getNodeReadCount(getTwinNode(node), graph);
	for (shortReadIndex = 0; shortReadIndex < shortReadCount; shortReadIndex++) {
		shortReadMarker = getShortReadMarkerAtIndex(shortReadArray, shortReadIndex);
		readID = getShortReadMarkerID(shortReadMarker);
		used[readID] = true;	
	}

}

void exportUnusedTranscriptReads(Graph* graph, Locus * loci, IDnum locusCount, ReadSet * reads, Coordinate minTransLength, char* directory) {
	char *outFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	FILE * outfile;
	boolean * used = callocOrExit(sequenceCount(graph) + 1, boolean);
	IDnum nodeID, readID;
	IDnum locusIndex;
	Transcript * transcript;

	strcpy(outFilename, directory);
	strcat(outFilename, "/UnusedReads.fa");
	outfile = fopen(outFilename, "w");

	velvetLog("Printing unused reads into %s\n", outFilename);

	resetNodeStatus(graph);

	for (locusIndex = 0; locusIndex < locusCount; locusIndex++)
		for (transcript = getTranscript(getLocus(loci, locusIndex)); transcript; transcript = getNextTranscript(transcript))
			if (getTranscriptLength(transcript) >= minTransLength)
				for(nodeID = 0; nodeID < getTranscriptContigCount(transcript); nodeID++)
					markUsedReads(getTranscriptContig(transcript, nodeID), used);

	for (readID = 1; readID <= sequenceCount(graph); readID++) 
		if (!used[readID])
			exportTightString(outfile, getTightStringInArray(reads->tSequences, readID - 1), readID);	

	free(outFilename);
	free(used);	
	fclose(outfile);
}

IDnum usedTranscriptReads(Graph * graph, Coordinate minTransLength, Locus * loci, IDnum locusCount) 
{
	boolean * used = callocOrExit(sequenceCount(graph) + 1, boolean);
	IDnum nodeID, readID;
	IDnum locusIndex;
	Transcript * transcript;
	IDnum res = 0;

	resetNodeStatus(graph);

	for (locusIndex = 0; locusIndex < locusCount; locusIndex++)
		for (transcript = getTranscript(getLocus(loci, locusIndex)); transcript; transcript = getNextTranscript(transcript))
			if (getTranscriptLength(transcript) >= minTransLength)
				for(nodeID = 0; nodeID < getTranscriptContigCount(transcript); nodeID++)
					markUsedReads(getTranscriptContig(transcript, nodeID), used);

	for (readID = 1; readID <= sequenceCount(graph); readID++) 
		if (used[readID])
			res++;

	free(used);	

	return res;
}

static IDnum getReferenceCount(ReadSet * reads) {
	IDnum index;

	for (index = 0; index < reads->readCount; index++) 
		if (reads->categories[index] <= 2 * CATEGORIES + 1)
			break;

	return index;
}

//////////////////////////////////////////////////////////////////////////
// Reference identifiers
//////////////////////////////////////////////////////////////////////////

typedef struct referenceCoord_st ReferenceCoord;

struct referenceCoord_st {
	char * name;
	Coordinate start;
	Coordinate finish;
	boolean positive_strand;
};

static ReferenceCoord * collectReferenceCoords(char * sequencesFilename, IDnum referenceCount) {
	FILE * file = fopen(sequencesFilename, "r");
	char line[MAXLINE];
	char name[500];
	Coordinate start, finish;
	long long longlongvar;
	IDnum refIndex = 0;
	ReferenceCoord * refCoords = callocOrExit(referenceCount, ReferenceCoord);
	int i;

	while (fgets(line, MAXLINE, file)) {
		if (line[0] == '>') {
			if (strchr(line, ':')) {
				sscanf(strtok(line, ":-\r\n"), ">%s", name);
				sscanf(strtok(NULL, ":-\r\n"), "%lli", &longlongvar);
				start = longlongvar;
				sscanf(strtok(NULL, ":-\r\n"), "%lli", &longlongvar);
				finish = longlongvar;
				refCoords[refIndex].name = callocOrExit(strlen(name) + 1, char);  
				if (start <= finish) {
					strcpy(refCoords[refIndex].name, name);
					refCoords[refIndex].start = start;
					refCoords[refIndex].finish = finish;
					refCoords[refIndex].positive_strand = true;
				} else {
					strcpy(refCoords[refIndex].name, name);
					refCoords[refIndex].start = finish;
					refCoords[refIndex].finish = start;
					refCoords[refIndex].positive_strand = false;
				}
			} else {
				for (i = strlen(line) - 1;
				     i >= 0 && (line[i] == '\n' || line[i] == '\r'); i--) {
					line[i] = '\0';
				}

				refCoords[refIndex].name = callocOrExit(strlen(name) + 1, char);  
				strcpy(name, line + 1);
				strcpy(refCoords[refIndex].name, name);
				refCoords[refIndex].start = 1;
				refCoords[refIndex].finish = -1;
				refCoords[refIndex].positive_strand = true;
			}
			if (++refIndex == referenceCount)
				break;	
		}
	}
	
	fclose(file);
	return refCoords;
}

typedef struct refMap_st {
	Coordinate start;
	Coordinate finish;
	IDnum refID;
	Coordinate refStart;
	Coordinate refFinish;
} ReferenceMapping; 

static int compareReferenceMappings(const void * A, const void * B) {
	ReferenceMapping * refMapA = (ReferenceMapping *) A;
	ReferenceMapping * refMapB = (ReferenceMapping *) B;
	
	if (refMapA->start < refMapB->start)
		return -1;
	else if (refMapA->start == refMapB->start)
		return 0;
	else 
		return 1;
}

static void initializeReferenceMapping(ReferenceMapping * refMap, PassageMarkerI marker, Transcript * transcript, IDnum nodeIndex, Coordinate nodeOffset) {
	PassageMarkerI finishMarker = marker;
	Coordinate totalLength = getNodeLength(getTranscriptContig(transcript, nodeIndex));
	IDnum index;

	for (index = nodeIndex + 1; index < getTranscriptContigCount(transcript); index++) {
		if (!getNextInSequence(finishMarker))
			break;

		if (getNode(getNextInSequence(finishMarker)) == getTranscriptContig(transcript, index)) {
			finishMarker = getNextInSequence(finishMarker);
			totalLength += getNodeLength(getTranscriptContig(transcript, index));
		}
	}

	refMap->start = nodeOffset + getStartOffset(marker);
	refMap->finish = nodeOffset + totalLength - getFinishOffset(finishMarker); 
	refMap->refID = getPassageMarkerSequenceID(marker);
	refMap->refStart = getPassageMarkerStart(marker);
	refMap->refFinish = getPassageMarkerFinish(finishMarker);
}

static void fprintfReferenceMapping(FILE * file, ReferenceMapping * mapping, ReferenceCoord * refCoords, int wordLength) {
	ReferenceCoord * refCoord;
	Coordinate start, finish;

	if (mapping->refID > 0) 
		refCoord = &refCoords[mapping->refID - 1];
	else
		refCoord = &refCoords[-mapping->refID - 1];

	if (mapping->refID > 0) {
		if (refCoord->positive_strand) {
			start = refCoord->start + mapping->refStart;
			finish = refCoord->start + mapping->refFinish + wordLength - 2;
		} else {
			start = refCoord->finish - mapping->refStart + wordLength - 1;
			finish = refCoord->finish - mapping->refFinish + 1;
		}
	} else {
		if (refCoord->positive_strand) {
			start = refCoord->start + mapping->refStart + wordLength - 1;
			finish = refCoord->start + mapping->refFinish + 1;
		} else {
			start = refCoord->finish - mapping->refStart; 
			finish = refCoord->finish - mapping->refFinish + wordLength;  
		}
	}
		
	fprintf(file, "%lli\t%lli\t%s\t%lli\t%lli\n",
		(long long) mapping->start + 1, (long long) mapping->finish + wordLength - 1, 
		refCoord->name, (long long) start, (long long) finish);
}

static void exportTranscriptMapping(FILE * outfile, Transcript * transcript, IDnum locusIndex, IDnum transcriptIndex, ReadSet * reads, ReferenceCoord * refCoords, int wordLength) {
	PassageMarkerI marker;
	ReferenceMapping * referenceMappings;
	IDnum index;
	IDnum referenceCount = 0;
	IDnum nodeIndex;
	Coordinate nodeOffset = 0;

	// Count reference sequences
	for (nodeIndex = 0; nodeIndex < getTranscriptContigCount(transcript); nodeIndex++)
		for (marker = getMarker(getTranscriptContig(transcript, nodeIndex)); marker; marker = getNextInNode(marker))
			if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] > 2 * CATEGORIES + 1
			    && (nodeIndex == 0 
				|| getNode(getPreviousInSequence(marker)) != getTranscriptContig(transcript, nodeIndex - 1)))
				referenceCount++;

	// Header
	fprintf(outfile, ">Locus_%li_Transcript_%li\n", (long) locusIndex + 1, (long) transcriptIndex + 1);

	// Create table
	referenceMappings = callocOrExit(referenceCount, ReferenceMapping);	

	// Initialize table
	referenceCount = 0;
	for (nodeIndex = 0; nodeIndex < getTranscriptContigCount(transcript); nodeIndex++) {
		for (marker = getMarker(getTranscriptContig(transcript, nodeIndex)); marker; marker = getNextInNode(marker))
			if (reads->categories[getAbsolutePassMarkerSeqID(marker) - 1] > 2 * CATEGORIES + 1
			    && (nodeIndex == 0 
				|| getNode(getPreviousInSequence(marker)) != getTranscriptContig(transcript, nodeIndex - 1)))
				initializeReferenceMapping(&referenceMappings[referenceCount++], marker, transcript, nodeIndex, nodeOffset);
		nodeOffset += getNodeLength(getTranscriptContig(transcript, nodeIndex));
		nodeOffset += getTranscriptDistance(transcript, nodeIndex);
	}

	// Sort table
	qsort(referenceMappings, referenceCount, sizeof(ReferenceMapping), compareReferenceMappings);

	// Print table
	for (index = 0; index < referenceCount; index++)
		fprintfReferenceMapping(outfile, &referenceMappings[index], refCoords, wordLength);

	// Clean table
	free(referenceMappings);
}

static void exportLocusMapping(FILE * outfile, Locus * loci, IDnum locusIndex, ReadSet * reads, ReferenceCoord * refCoords, Coordinate minTransLength, int wordLength) {
	Transcript * transcript;
	IDnum transcriptIndex = 0;

	for (transcript = getTranscript(getLocus(loci, locusIndex)); transcript != NULL;
	     transcript = getNextTranscript(transcript))
		if (getTranscriptLength(transcript) >= minTransLength)
			exportTranscriptMapping(outfile, transcript, locusIndex, transcriptIndex++, reads, refCoords, wordLength);
		
}

void exportTranscriptMappings(Locus * loci, IDnum locusCount, 
			      Graph * graph, ReadSet * reads,
			      Coordinate minLength, char * directory)
{
	FILE * outfile;
	IDnum locusIndex, refIndex;
	ReferenceCoord * refCoords;
	IDnum referenceCount = getReferenceCount(reads); 
	char filename[10000];

	strcpy(filename, directory);
	strcat(filename, "/transcript-alignments.psa");
	velvetLog("Writing into pseudo-alignment file %s...\n", filename);
	outfile = fopen(filename, "w");

	if (referenceCount == 0)	
		return;

	strcpy(filename, directory);
	strcat(filename, "/Sequences");
	refCoords = collectReferenceCoords(filename, referenceCount);

	strcpy(filename, directory);
	strcat(filename, "/transcript-alignments.psa");
	outfile = fopen(filename, "w");
	if (outfile == NULL) {
		velvetLog("Could not write into %s, sorry\n", filename);
		return;
	} else {
		velvetLog("Writing contigs into %s...\n", filename);
	}

	for (locusIndex = 0; locusIndex < locusCount; locusIndex++) 
		exportLocusMapping(outfile, loci, locusIndex, reads, refCoords, minLength, getWordLength(graph));

	for (refIndex = 0; refIndex < referenceCount; refIndex++)
		free(refCoords[refIndex].name);
	free(refCoords);
	fclose(outfile);
}


static void exportTranscript(Transcript * transcript, IDnum locusID,
			     IDnum transID, IDnum transcriptCount, FILE * outfile)
{
	IDnum index;
	Coordinate offset = getWordLength(graph) - 1;
	char * sequence = nSequence(getTranscriptLength(transcript)); 

	// Extract sequence
	for (index = 0; index < getTranscriptContigCount(transcript); index++) {
	    addNodeToSequence(sequence, getTranscriptContig(transcript, index), offset);

	    // Increment offset
	    offset += getNodeLength(getTranscriptContig(transcript, index));
	    if (index < getTranscriptContigCount(transcript) - 1)
		    offset += getTranscriptDistance(transcript, index);
	}

	// Count initial N's
	for (offset = 0; offset < strlen(sequence); offset++)
		if (sequence[offset] != 'N')
			break;
	
	// Print
	fprintf(outfile, ">Locus_%li_Transcript_%li/%li_Confidence_%.3f_Length_%li\n",
		(long) locusID + 1, (long) transID + 1, (long) transcriptCount, getConfidence(transcript), (long) (strlen(sequence) - offset));
	printFastA(outfile, sequence + offset);

	free(sequence);
}

static void exportLocusTranscripts(Locus * locus, IDnum locusID,
				   FILE * outfile, Coordinate minTransLength)
{
	IDnum index = 0;
	Transcript *transcript;
	IDnum transcriptCount = 0;

	for (transcript = getTranscript(locus); transcript != NULL;
	     transcript = getNextTranscript(transcript))
		if (getTranscriptLength(transcript) >= minTransLength)
			transcriptCount++;

	for (transcript = getTranscript(locus); transcript != NULL;
	     transcript = getNextTranscript(transcript))
		if (getTranscriptLength(transcript) >= minTransLength)
			exportTranscript(transcript, locusID, index++, transcriptCount, outfile);
}

void exportTranscripts(Locus * loci, IDnum locusCount, char *filename, Coordinate minTransLength, Graph * argGraph)
{
	FILE *outfile = fopen(filename, "w");
	IDnum index;
	graph = argGraph;

	velvetLog("Exporting transcripts to %s\n", filename);

	for (index = 0; index < locusCount; index++)
		exportLocusTranscripts(getLocus(loci, index), index, outfile, minTransLength);
	fclose(outfile);
}

void logFinalOasesStats(Graph * graph, Coordinate minTransLength, Locus * loci, IDnum locusCount, char *directory) {
	char *logFilename =
	    mallocOrExit(strlen(directory) + 100, char);
	char *statsLine = 
	    mallocOrExit(5000, char);
	FILE *logFile;

	strcpy(logFilename, directory);
	strcat(logFilename, "/Log");
	logFile = fopen(logFilename, "a");

	if (logFile == NULL)
		exitErrorf(EXIT_FAILURE, true, "Could not write to %s",
		       logFilename);

	sprintf
	    (statsLine, "Finished extracting transcripts, used %li/%li reads on %li loci\n", (long) usedTranscriptReads(graph, minTransLength, loci, locusCount), (long) sequenceCount(graph), (long) locusCount);

	velvetFprintf(logFile, "%s", statsLine);
	velvetLog("%s", statsLine);

	fclose(logFile);
	free(logFilename);
	free(statsLine);
}
