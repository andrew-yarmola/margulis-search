#REFINEOBJS= Refine.o TubeSearch.o TestCollection.o Box.o ImpossibleRelations.o AJ.o SL2C.o SL2AJ.o Complex.o types.o roundoff.o CanonicalName.o QuasiRelators.o
REFINEOBJS= Refine.o Box.o AJ.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o ImpossibleRelations.o QuasiRelators.o
TESTOBJS= test.o Box.o AJ.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o
CENSUSTESTOBJS= census_test.o Box.o AJ.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o ImpossibleRelations.o QuasiRelators.o
PLAYOBJS= playground.o Box.o AJ.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o
TREECATOBJS= treecat.o
CPLUSPLUS= g++
CC= g++

#COPTS=-g3 -std=c++11 -stdlib=libc++
COPTS=-O3 -std=c++1y

.cpp.o:
	$(CPLUSPLUS) $(COPTS) $(CCFLAGS) -c ${@:o=cpp}

.C.o:
	$(CC) $(COPTS) $(CCFLAGS) -c ${@:o=C}

test: $(TESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(TESTOBJS) $(LIBS) -lm -o test

censustest: $(CENSUSTESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(CENSUSTESTOBJS) $(LIBS) -lm -o census_test

playground: $(TESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(PLAYOBJS) $(LIBS) -lm -o playground

refine: $(REFINEOBJS)
	$(CPLUSPLUS) $(COPTS) $(REFINEOBJS) $(LIBS) -lm -o ../bin/refine_marg

treecat: $(TREECATOBJS)
	$(CPLUSPLUS) $(COPTS) $(TREECATOBJS) -o ../bin/treecat
