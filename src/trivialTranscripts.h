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
#ifndef _TRIVIALTRANSCRIPTS_H_
#define _TRIVIALTRANSCRIPTS_H_

#define UNKNOWN 0
#define LINEAR 1
#define FORK 2
#define BUBBLE 3

void getDegreeDistribution(Locus * locus, IDnum * distribution);
boolean hasPlausibleTranscripts(IDnum * distribution);
void addLinearTranscriptToLocus(Locus * locus);
void addBubbleTranscriptsToLocus(Locus * locus);
void addForkedTranscriptsToLocus(Locus * locus);
#endif
