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
#ifndef _SCAFFOLD_H_
#define _SCAFFOLD_H_

typedef struct connection_st Connection;

void buildScaffold(Graph * graph, ReadSet * reads, boolean * dubious,
		   IDnum * lengths, boolean scaffolding);
void transitiveReduction();

Connection *getConnection(Node * node);
Connection *getConnectionBetweenNodes(Node * nodeA, Node * nodeB);
Connection *getNextConnection(Connection * connect);

boolean getConnectionStatus(Connection * connect);
void setConnectionStatus(Connection * connect, boolean status);
Node *getConnectionDestination(Connection * connect);
void incrementConnectionWeight(Connection * connect, double weight);
double getConnectionWeight(Connection * connect);
void setConnectionWeight(Connection * connect, double weight);
Coordinate getConnectionDistance(Connection * connect);
IDnum getConnectionDirectCount(Connection * connect);

void printOasesConnections(Category * cats);
void scaffold_setPairedThreshold(double pairedThreshold);
void scaffold_setUnreliableConnectionCutoff(int val);
void scaffold_setDegreeCutoff(int val);

void destroyConnection(Connection * connect, IDnum index);
void cleanScaffoldMemory();
#endif
