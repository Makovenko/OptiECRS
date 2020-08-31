
OBJS	= optimizemixedinteger.o optimizebranchandbound.o optimizedirection.o randomheuristic.o optimizeconsecutive.o optimizedirectionalcost.o optimizemixedintegeralt.o galois.o gf_general.o gf_cpu.o gf_w4.o cauchy.o gf.o liberation.o gf_w32.o gf_wgen.o gf_w8.o galois_d.o gf_w16.o gf_method.o gf_rand.o timing.o jerasure.o gf_w64.o reed_sol.o gf_w128.o extendedcauchymatrix.o main.o
SOURCE	= procedures/optimizemixedinteger.cpp procedures/optimizebranchandbound.cpp procedures/optimizedirection.cpp procedures/randomheuristic.cpp procedures/optimizeconsecutive.cpp procedures/optimizedirectionalcost.cpp procedures/optimizemixedintegeralt.cpp galois/galois.cpp external/bundle/gf_general.c external/bundle/gf_cpu.c external/bundle/gf_w4.c external/bundle/cauchy.c external/bundle/gf.c external/bundle/liberation.c external/bundle/gf_w32.c external/bundle/gf_wgen.c external/bundle/gf_w8.c external/bundle/galois.c external/bundle/gf_w16.c external/bundle/gf_method.c external/bundle/gf_rand.c external/bundle/timing.c external/bundle/jerasure.c external/bundle/gf_w64.c external/bundle/reed_sol.c external/bundle/gf_w128.c matrix/extendedcauchymatrix.cpp main.cpp
HEADER	= procedures/optimizedirection.h procedures/optimizemixedinteger.h procedures/randomheuristic.h procedures/optimizemixedintegeralt.h procedures/optimizebranchandbound.h procedures/optimizedirectionalcost.h procedures/optimizeconsecutive.h galois/galois.h external/bundle/gf_w16.h external/bundle/gf_method.h external/bundle/galois.h external/bundle/gf_complete.h external/bundle/gf_w8.h external/bundle/gf_rand.h external/bundle/jerasure.h external/bundle/gf_int.h external/bundle/timing.h external/bundle/reed_sol.h external/bundle/gf_w64.h external/bundle/gf_w4.h external/bundle/gf_general.h external/bundle/gf_cpu.h external/bundle/cauchy.h external/bundle/liberation.h external/bundle/gf_w32.h matrix/extendedcauchymatrix.h
OUT	= OptiECRS
CC	 = g++
C	 = gcc
FLAGS	 = -c -Wall -O3
GUROBIINC = -I/Library/gurobi903/mac64/include
GUROBILIB = -L/Library/gurobi903/mac64/lib
LFLAGS	 = 

all: $(OBJS)
	$(CC) -g $(OBJS) -o $(OUT) $(LFLAGS) $(GUROBIINC) $(GUROBILIB) -lgurobi90 -lgurobi_c++

optimizemixedinteger.o: procedures/optimizemixedinteger.cpp
	$(CC) $(FLAGS) $(GUROBIINC) procedures/optimizemixedinteger.cpp -std=c++17

optimizebranchandbound.o: procedures/optimizebranchandbound.cpp
	$(CC) $(FLAGS) procedures/optimizebranchandbound.cpp -std=c++17

optimizedirection.o: procedures/optimizedirection.cpp
	$(CC) $(FLAGS) procedures/optimizedirection.cpp -std=c++17

randomheuristic.o: procedures/randomheuristic.cpp
	$(CC) $(FLAGS) procedures/randomheuristic.cpp -std=c++17

optimizeconsecutive.o: procedures/optimizeconsecutive.cpp
	$(CC) $(FLAGS) procedures/optimizeconsecutive.cpp -std=c++17

optimizedirectionalcost.o: procedures/optimizedirectionalcost.cpp
	$(CC) $(FLAGS) procedures/optimizedirectionalcost.cpp -std=c++17

optimizemixedintegeralt.o: procedures/optimizemixedintegeralt.cpp
	$(CC) $(FLAGS) $(GUROBIINC) procedures/optimizemixedintegeralt.cpp -std=c++17

galois_d.o: galois/galois.cpp
	$(CC) $(FLAGS) galois/galois.cpp -std=c++17 -o galois_d.o

gf_general.o: external/bundle/gf_general.c
	$(C) $(FLAGS) external/bundle/gf_general.c -std=c17

gf_cpu.o: external/bundle/gf_cpu.c
	$(C) $(FLAGS) external/bundle/gf_cpu.c -std=c17

gf_w4.o: external/bundle/gf_w4.c
	$(C) $(FLAGS) external/bundle/gf_w4.c -std=c17

cauchy.o: external/bundle/cauchy.c
	$(C) $(FLAGS) external/bundle/cauchy.c -std=c17

gf.o: external/bundle/gf.c
	$(C) $(FLAGS) external/bundle/gf.c -std=c17

liberation.o: external/bundle/liberation.c
	$(C) $(FLAGS) external/bundle/liberation.c -std=c17

gf_w32.o: external/bundle/gf_w32.c
	$(C) $(FLAGS) external/bundle/gf_w32.c -std=c17

gf_wgen.o: external/bundle/gf_wgen.c
	$(C) $(FLAGS) external/bundle/gf_wgen.c -std=c17

gf_w8.o: external/bundle/gf_w8.c
	$(C) $(FLAGS) external/bundle/gf_w8.c -std=c17

galois.o: external/bundle/galois.c
	$(C) $(FLAGS) external/bundle/galois.c -std=c17

gf_w16.o: external/bundle/gf_w16.c
	$(C) $(FLAGS) external/bundle/gf_w16.c -std=c17

gf_method.o: external/bundle/gf_method.c
	$(C) $(FLAGS) external/bundle/gf_method.c -std=c17

gf_rand.o: external/bundle/gf_rand.c
	$(C) $(FLAGS) external/bundle/gf_rand.c -std=c17
r6.c -std=c17

timing.o: external/bundle/timing.c
	$(C) $(FLAGS) external/bundle/timing.c -std=c17

jerasure.o: external/bundle/jerasure.c
	$(C) $(FLAGS) external/bundle/jerasure.c -std=c17

gf_w64.o: external/bundle/gf_w64.c
	$(C) $(FLAGS) external/bundle/gf_w64.c -std=c17

reed_sol.o: external/bundle/reed_sol.c
	$(C) $(FLAGS) external/bundle/reed_sol.c -std=c17

gf_w128.o: external/bundle/gf_w128.c
	$(C) $(FLAGS) external/bundle/gf_w128.c -std=c17

extendedcauchymatrix.o: matrix/extendedcauchymatrix.cpp
	$(CC) $(FLAGS) matrix/extendedcauchymatrix.cpp -std=c++17

main.o: main.cpp
	$(CC) $(FLAGS) $(GUROBIINC) main.cpp -std=c++17


clean:
	rm -f $(OBJS) $(OUT)
	rm -f $(OUT)
