include Makefile.inc
.DEFAULT_GOAL = randomice
PROG = randomice

CCFLGS = 
LDFLAGS = 
HEADERS = atom.h cell.h randomice.h
SRCS = atom.cc cell.cc randomice.cc
OBJS = atom.o cell.o randomice.o

.SUFFIXES: $(SUFFIXES) .cc .h
.PHONY: clean

all: $(PROG)

clean:
	rm -f $(PROG) $(OBJS)

$(PROG): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

.cc.o:
	$(CC) -c $< $(LDFLAGS)

atom.o: atom.h atom.cc
cell.o: atom.h atom.cc cell.h cell.cc
randomice.o: atom.h atom.cc cell.h cell.cc randomice.h randomice.cc
