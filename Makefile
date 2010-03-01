CC = gcc
CFLAGS = -Wall
DEBUG = -g
LDFLAGS = -lm
OPT = -O3
MAXKMERLENGTH=31
CATEGORIES=2
DEF = -D MAXKMERLENGTH=$(MAXKMERLENGTH) -D CATEGORIES=$(CATEGORIES)

VELVET_DIR=../../velvet
VELVET_SRC_DIR=$(VELVET_DIR)/src
VELVET_OBJ = recycleBin utility graph passageMarker readSet tightString kmer dfibHeap dfib concatenatedGraph 
VELVET_FILES = $(VELVET_OBJ:%=$(VELVET_DIR)/obj/%.o)
VELVET_DBG_FILES = $(VELVET_OBJ:%=$(VELVET_DIR)/obj/dbg/%.o)

Z_LIB_DIR=$(VELVET_DIR)/third-party/zlib-1.2.3
Z_LIB_FILES=$(Z_LIB_DIR)/*.o

# Mac OS users: uncomment the following lines
# Z_LIB_FILES=
# LDFLAGS = -lm -lz
# CFLAGS = -Wall -m64

# Sparc/Solaris users: uncomment the following line
# CFLAGS = -Wall -m64

OBJ = obj/oases.o obj/transcript.o obj/scaffold.o
OBJDBG = $(subst obj,obj/dbg,$(OBJ))

default : oases

velvet :
	cd $(VELVET_DIR) && make obj

velvetdbg :
	cd $(VELVET_DIR) && make obj/dbg

velvet_de :
	cd $(VELVET_DIR) && make obj_de

velvetdbg_de :
	cd $(VELVET_DIR) && make obj/dbg_de

clean :
	-rm obj/*.o obj/dbg/*.o ./oases 
	cd $(VELVET_DIR) && make clean

cleanobj: 
	-rm obj/*.o obj/dbg/*.o 

oases : cleanobj velvet obj $(OBJ) 
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o oases $(OBJ) $(VELVET_FILES) $(Z_LIB_FILES)


debug : cleanobj velvetdbg obj/dbg $(OBJDBG)
	$(CC) $(CFLAGS) $(DEBUG) $(LDFLAGS) -o oases $(OBJDBG) $(VELVET_DBG_FILES) $(Z_LIB_FILES)

color : override DEF := $(DEF) -D COLOR
color : cleanobj velvet_de obj $(OBJ)
	$(CC) $(CFLAGS) $(OPT) $(LDFLAGS) -o oases_de $(OBJ) $(VELVET_FILES) $(Z_LIB_FILES)

colordebug : override DEF := $(DEF) -D COLOR
colordebug : cleanobj velvetdbg_de obj/dbg $(OBJDBG) 
	$(CC) $(CFLAGS) $(DEBUG) $(LDFLAGS) -o oases_de $(OBJDBG) $(VELVET_DBG_FILES) $(Z_LIB_FILES)

obj:
	mkdir -p obj

obj/dbg: 
	mkdir -p obj/dbg

obj/%.o: src/%.c
	$(CC) $(CFLAGS) $(OPT) $(DEF) -c $? -o $@ -I$(VELVET_SRC_DIR)

obj/dbg/%.o: src/%.c
	$(CC) $(CFLAGS) $(DEBUG) $(DEF) -c $? -o $@ -I$(VELVET_SRC_DIR)
