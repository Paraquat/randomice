include Makefile.inc
.DEFAULT_GOAL = randomice
PROG = randomice

LIBS = -lboost_program_options
CCFLGS = 
LDFLAGS = 
HEADERS = atom.h cell.h ice.h randomice.h
SRCS = atom.cc cell.cc ice.cc randomice.cc
OBJS = atom.o cell.o ice.o randomice.o

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
cell.o: atom.h atom.cc cell.h cell.cc
ice.o: atom.h atom.cc cell.h cell.cc ice.h ice.cc
randomice.o: atom.h atom.cc cell.h cell.cc ice.h ice.cc randomice.h randomice.cc