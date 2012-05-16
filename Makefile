CC = gcc
CFLAGS = -Wall
DEBUG = -g
LIBS = -lm
OPT = -O3
export MAXKMERLENGTH = 64
export CATEGORIES = 2
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH) -D CATEGORIES=$(CATEGORIES)
VELVET_DIR=../velvet

ifdef BIGASSEMBLY
override DEF := $(DEF) -D BIGASSEMBLY
endif 	

ifdef VBIGASSEMBLY
override DEF := $(DEF) -D BIGASSEMBLY -D VBIGASSEMBLY
endif 	

ifdef LONGSEQUENCES
override DEF := $(DEF) -D LONGSEQUENCES
endif 	

ifdef OPENMP
override CFLAGS := $(CFLAGS) -fopenmp
endif

# Per library coverage
ifdef SINGLE_COV_CAT
override DEF := $(DEF) -D SINGLE_COV_CAT
endif

default : cleanobj velvet oases doc

ifdef BUNDLEDZLIB
Z_LIB_DIR=$(VELVET_DIR)/third-party/zlib-1.2.3
Z_LIB_FILES=$(Z_LIB_DIR)/*.o
override DEF := $(DEF) -D BUNDLEDZLIB

zlib: 
	cd $(Z_LIB_DIR); ./configure; make; rm minigzip.o; rm example.o

clean-zlib:
	cd $(Z_LIB_DIR) && make clean
else
Z_LIB_FILES=-lz

zlib:
clean-zlib:
endif

VELVET_SRC_DIR=$(VELVET_DIR)/src
VELVET_OBJ = recycleBin utility graph passageMarker readSet tightString kmer dfibHeap dfib concatenatedGraph graphStats fibHeap fib readCoherentGraph allocArray binarySequences
VELVET_FILES = $(VELVET_OBJ:%=$(VELVET_DIR)/obj/%.o)
VELVET_DBG_FILES = $(VELVET_OBJ:%=$(VELVET_DIR)/obj/dbg/%.o)

# Mac OS users: uncomment the following lines
# Z_LIB_FILES=
# LDFLAGS = -lm -lz
# CFLAGS = -Wall -m64

# Sparc/Solaris users: uncomment the following line
# CFLAGS = -Wall -m64

OBJ = obj/oases.o obj/transcript.o obj/scaffold.o obj/locallyCorrectedGraph2.o obj/correctedGraph.o obj/filterTranscripts.o obj/locus.o obj/nodeList.o obj/oasesExport.o obj/trivialTranscripts.o obj/complexTranscripts.o obj/extractMergedTranscripts.o obj/extractLoci.o
OBJDBG = $(subst obj,obj/dbg,$(OBJ))

velvet :
	cd $(VELVET_DIR) && make -e obj

velvetdbg :
	cd $(VELVET_DIR) && make -e obj/dbg

velvet_de :
	cd $(VELVET_DIR) && make -e obj_de

velvetdbg_de :
	cd $(VELVET_DIR) && make -e obj/dbg_de

clean :
	rm -f obj/*.o obj/dbg/*.o ./oases 
	cd $(VELVET_DIR) && make clean
	cd doc && make clean

cleanobj: 
	rm -f obj/*.o obj/dbg/*.o 

doc: OasesManual.pdf

OasesManual.pdf: doc/manual/OasesManual.tex
	cd doc; make

oases : obj $(OBJ) 
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o oases $(OBJ) $(VELVET_FILES) $(Z_LIB_FILES) $(LIBS)


debug : cleanobj velvetdbg obj/dbg $(OBJDBG)
	$(CC) $(CFLAGS) $(DEBUG) $(LDFLAGS) -o oases $(OBJDBG) $(VELVET_DBG_FILES) $(Z_LIB_FILES) $(LIBS)

color : override DEF := $(DEF) -D COLOR
color : cleanobj velvet_de obj $(OBJ)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o oases_de $(OBJ) $(VELVET_FILES) $(Z_LIB_FILES) $(LIBS)

colordebug : override DEF := $(DEF) -D COLOR
colordebug : cleanobj velvetdbg_de obj/dbg $(OBJDBG) 
	$(CC) $(CFLAGS) $(DEBUG) $(LDFLAGS) -o oases_de $(OBJDBG) $(VELVET_DBG_FILES) $(Z_LIB_FILES) $(LIBS)

obj:
	mkdir -p obj

obj/dbg: 
	mkdir -p obj/dbg

obj/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPT) $(DEF) -c $? -o $@ -I$(VELVET_SRC_DIR)

obj/dbg/%.o: src/%.c
	$(CC) $(CFLAGS) $(DEBUG) $(DEF) -c $? -o $@ -I$(VELVET_SRC_DIR)
