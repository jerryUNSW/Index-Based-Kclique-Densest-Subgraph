
OBJECT_DIR = obj
SRC_DIR    = src
BIN_DIR    = bin

OBJECTS  = $(OBJECT_DIR)/LinkedList.o
OBJECTS += $(OBJECT_DIR)/MemoryManager.o
OBJECTS += $(OBJECT_DIR)/degeneracy_algorithm_cliques_A.o
OBJECTS += $(OBJECT_DIR)/degeneracy_helper.o
OBJECTS += $(OBJECT_DIR)/misc.o

EXEC_NAMES =  compdegen degeneracy_cliques

EXECS = $(addprefix $(BIN_DIR)/, $(EXEC_NAMES))

VPATH = src

.PHONY : all

all: $(EXECS)

.PHONY : clean

clean: 
	rm -rf $(OBJECTS) $(EXECS) $(OBJECT_DIR) $(BIN_DIR)

$(BIN_DIR)/compdegen: compdegen.c ${OBJECTS} ${BIN_DIR}
	g++ -O3 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/compdegen.c -o $@ -fopenmp


$(BIN_DIR)/degeneracy_cliques: degeneracy_cliques.c ${OBJECTS} ${BIN_DIR}
	g++ -O3 -g ${DEFINE} ${OBJECTS} $(SRC_DIR)/degeneracy_cliques.c -o $@ -fopenmp

$(OBJECT_DIR)/LinkedList.o: LinkedList.c LinkedList.h ${OBJECT_DIR}
	g++ -O3 -g ${DEFINE} -c $(SRC_DIR)/LinkedList.c -o $@ -fopenmp

$(OBJECT_DIR)/MemoryManager.o: MemoryManager.c MemoryManager.h ${OBJECT_DIR}
	g++ -O3 -g ${DEFINE} -c $(SRC_DIR)/MemoryManager.c -o $@ -fopenmp

$(OBJECT_DIR)/degeneracy_algorithm_cliques_A.o: degeneracy_algorithm_cliques_A.c degeneracy_algorithm_cliques_A.h ${OBJECT_DIR}
	g++ -O3 -g ${DEFINE} -c $(SRC_DIR)/degeneracy_algorithm_cliques_A.c -o $@ -fopenmp

$(OBJECT_DIR)/degeneracy_helper.o: degeneracy_helper.c degeneracy_helper.h ${OBJECT_DIR}
	g++ -O3 -g ${DEFINE} -c $(SRC_DIR)/degeneracy_helper.c -o $@ -fopenmp

$(OBJECT_DIR)/misc.o: misc.c misc.h ${OBJECT_DIR}
	g++ -O3 -g ${DEFINE} -c $(SRC_DIR)/misc.c -o $@ -fopenmp

${OBJECT_DIR}:
	mkdir ${OBJECT_DIR}

${BIN_DIR}:
	mkdir ${BIN_DIR}
