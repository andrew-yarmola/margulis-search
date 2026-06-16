REFINEOBJS= RefineMain.o Refine.o PartialTree.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o QuasiRelators.o RelatorTest.o
SEARCHTESTOBJS= TestMain.o Testing.o PartialTree.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o QuasiRelators.o RelatorTest.o
TESTOBJS= test.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o
CENSUSTESTOBJS= census_test.o Refine.o PartialTree.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o TestCollection.o QuasiRelators.o RelatorTest.o
PLAYOBJS= playground.o Box.o AJCC.o Complex.o roundoff.o types.o TubeSearch.o CanonicalName.o
RECEDETESTOBJS= test_recede.o Box.o AJCC.o Complex.o roundoff.o types.o CanonicalName.o QuasiRelators.o
TREECATOBJS= treecat.o
TREESTATOBJS= treestat.o
CPLUSPLUS= g++ -std=c++1y
CC= g++ -std=c++1y

#COPTS=-g -Wall
# -march=native: tune for the build CPU. On a heterogeneous cluster, pin a specific
#   baseline (e.g. -march=x86-64-v3) so binaries run on every node.
# -ffp-contract=off: forbid FMA contraction so results are bit-reproducible across
#   machines/compilers, which a verified computation needs. Never add -ffast-math:
#   it reorders FP ops and breaks the AJCC (1+EPS) error model.
COPTS=-O3 -march=native -ffp-contract=off


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

# Soundness regression test for the subtree word-pruning predicate
# (word_provably_receding). Run: ../bin/recede_test   (sweeps deep boxes, exit 0 = clean)
recede_test: $(RECEDETESTOBJS)
	$(CPLUSPLUS) $(COPTS) $(RECEDETESTOBJS) $(LIBS) -lm -o ../bin/recede_test

refine: $(REFINEOBJS)
	$(CPLUSPLUS) $(COPTS) $(REFINEOBJS) $(LIBS) -lm -o ../bin/refine_marg

debug: $(REFINEOBJS)
	$(CPLUSPLUS) $(COPTS) $(REFINEOBJS) $(LIBS) -lm -o ../bin/refine_marg_debug

full: $(REFINEOBJS)
	$(CPLUSPLUS) $(COPTS) $(REFINEOBJS) $(LIBS) -lm -o ../bin/refine_marg_full

treecat: $(TREECATOBJS)
	$(CPLUSPLUS) $(COPTS) $(TREECATOBJS) -o ../bin/treecat

treestat: $(TREESTATOBJS)
	$(CPLUSPLUS) $(COPTS) $(TREESTATOBJS) -o ../bin/treestat

treemerge: treemerge.o 
	$(CPLUSPLUS) $(COPTS) treemerge.o -o ../bin/treemerge

roottreemerge: roottreemerge.o 
	$(CPLUSPLUS) $(COPTS) roottreemerge.o -o ../bin/roottreemerge
