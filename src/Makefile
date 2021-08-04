REFINEOBJS= RefineMain.o Refine.o PartialTree.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o QuasiRelators.o ImpossibleRelations.o
SEARCHTESTOBJS= TestMain.o Testing.o PartialTree.o Box.o AJ.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o QuasiRelators.o ImpossibleRelations.o
TESTOBJS= test.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o
CENSUSTESTOBJS= census_test.o Refine.o PartialTree.o Box.o AJ.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o QuasiRelators.o ImpossibleRelations.o
PLAYOBJS= playground.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o
TREECATOBJS= treecat.o
TREESTATOBJS= treestat.o
CPLUSPLUS= g++
CC= g++

#COPTS=-pg -Wall -std=c++1y
COPTS=-O3 -std=c++1y

.cpp.o:
	$(CPLUSPLUS) $(COPTS) $(CCFLAGS) -c ${@:o=cpp}

.C.o:
	$(CC) $(COPTS) $(CCFLAGS) -c ${@:o=C}

test: $(TESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(TESTOBJS) $(LIBS) -lm -o test

searchtest: $(SEARCHTESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(SEARCHTESTOBJS) $(LIBS) -lm -o ../bin/search_test

censustest: $(CENSUSTESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(CENSUSTESTOBJS) $(LIBS) -lm -o census_test

playground: $(PLAYOBJS)
	$(CPLUSPLUS) $(COPTS) $(PLAYOBJS) $(LIBS) -lm -o playground

refine: $(REFINEOBJS)
	$(CPLUSPLUS) $(COPTS) $(REFINEOBJS) $(LIBS) -lm -o ../bin/refine_marg

debug: $(REFINEOBJS)
	$(CPLUSPLUS) $(COPTS) $(REFINEOBJS) $(LIBS) -lm -o ../bin/refine_marg_debug

treecat: $(TREECATOBJS)
	$(CPLUSPLUS) $(COPTS) $(TREECATOBJS) -o ../bin/treecat

treestat: $(TREESTATOBJS)
	$(CPLUSPLUS) $(COPTS) $(TREESTATOBJS) -o ../bin/treestat
