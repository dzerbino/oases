CC = gcc
CFLAGS = -Wall
DEBUG = -g
LDFLAGS = -lm -lz
OPT = -O3
MAXKMERLENGTH=31
CATEGORIES=2
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH) -D CATEGORIES=$(CATEGORIES)

VELVET_DIR=
VELVET_SRC_DIR=$(VELVET_DIR)/src
VELVET_OBJ = recycleBin utility graph passageMarker readSet tightString kmer dfibHeap dfib concatenatedGraph 
VELVET_FILES = $(VELVET_OBJ:%=$(VELVET_DIR)/obj/%.o)
VELVET_DBG_FILES = $(VELVET_OBJ:%=$(VELVET_DIR)/obj/dbg/%.o)

# Mac OS users: uncomment the following lines
# Sparc/Solaris users: uncomment the following line
# CFLAGS = -Wall -m64

OBJ = obj/oases.o obj/transcript.o obj/scaffold.o
OBJDBG = $(subst obj,obj/dbg,$(OBJ))

default : cleanobj obj velvet oases

velvet :
	cd $(VELVET_DIR) && make

velvetdbg :
	cd $(VELVET_DIR) && make debug

clean :
	-rm obj/*.o obj/dbg/*.o ./oases 
	-rm -R obj/

cleanobj: 
	-rm obj/*.o obj/dbg/*.o 

oases : cleanobj obj $(OBJ) 
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o oases $(OBJ) $(VELVET_FILES)


debug : cleanobj velvetdbg obj/dbg $(OBJDBG)
	$(CC) $(CFLAGS) $(DEBUG) $(LDFLAGS) -o oases $(OBJDBG) $(VELVET_DBG_FILES)

color : override DEF := $(DEF) -D COLOR
color : cleanobj obj $(OBJ)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o oases_de $(OBJ) $(VELVET_FILES)

colordebug : override DEF := $(DEF) -D COLOR
colordebug : cleanobj velvetdbg obj/dbg $(OBJDBG) 
	$(CC) $(CFLAGS) $(DEBUG) $(LDFLAGS) -o oases_de $(OBJDBG) $(VELVET_DBG_FILES)

obj:
	mkdir -p obj

obj/dbg: 
	mkdir -p obj/dbg

obj/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPT) $(DEF) -c $? -o $@ -I$(VELVET_SRC_DIR)

obj/dbg/%.o: src/%.c
	$(CC) $(CFLAGS) $(DEBUG) $(DEF) -c $? -o $@ -I$(VELVET_SRC_DIR)
