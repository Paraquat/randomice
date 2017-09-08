include Makefile.inc
.DEFAULT_GOAL = randomice
PROG = randomice

LIBS = -lboost_program_options
CCFLGS = 
LDFLAGS = 
HEADERS = atom.h water.h hbond.h cell.h ice.h randomice.h constants.h
SRCS = atom.cc water.cc hbond.cc cell.cc ice.cc randomice.cc
OBJS = atom.o water.o hbond.o cell.o ice.o randomice.o

.SUFFIXES: $(SUFFIXES) .cc .h
.PHONY: clean

all: $(PROG)

clean:
	rm -f $(PROG) $(OBJS)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

.cc.o:
	$(CC) -c $< $(LDFLAGS)

atom.o: atom.h atom.cc
water.o: water.h water.cc atom.o
hbond.o: hbond.h hbond.cc atom.o
cell.o: atom.o cell.h cell.cc constants.h
ice.o: atom.o water.o hbond.o cell.o ice.h ice.cc
randomice.o: atom.o water.o hbond.o cell.o ice.o randomice.h randomice.cc