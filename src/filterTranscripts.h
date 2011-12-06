#ifndef FILTERTRANSCRIPTS_H_
#define FILTERTRASNCRIPTS_H_

void removeRedundantTranscripts(Graph * graph);
Locus *reextractGraphLoci(Graph * graph, IDnum * locusCount);
void recomputeTranscripts(Locus ** loci, IDnum * locusCount);
#endif 
