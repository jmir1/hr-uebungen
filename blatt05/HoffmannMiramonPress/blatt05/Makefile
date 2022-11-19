# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
LIBS   = -lm

OBJS = partdiff.o askparams.o
TGTS = partdiff-seq partdiff-openmp partdiff-openmp-zeile partdiff-openmp-spalte partdiff-openmp-element



# Rule to create *.o from *.c
.c.o:
	$(CC) -c $(CFLAGS) $*.c

# Targets ...
all: $(TGTS)
partdiff-seq: $(OBJS) Makefile
	$(CC) $(LFLAGS) -o $@ $(OBJS) $(LIBS)

partdiff-openmp: partdiff-openmp.o askparams.o Makefile
	gcc $(LFLAGS) -fopenmp -o $@ partdiff-openmp.o askparams.o -lm

partdiff-openmp-element: partdiff-openmp-element.o askparams.o Makefile
	gcc $(LFLAGS) -fopenmp -D ELEMENT -o $@ partdiff-openmp-element.o askparams.o -lm

partdiff-openmp-spalte: partdiff-openmp-spalte.o askparams.o Makefile
	gcc $(LFLAGS) -fopenmp -D SPALTE -o $@ partdiff-openmp-spalte.o askparams.o -lm

partdiff-openmp-zeile: partdiff-openmp-zeile.o askparams.o Makefile
	gcc $(LFLAGS) -fopenmp -D ZEILE -o $@ partdiff-openmp-zeile.o askparams.o -lm



partdiff.o: partdiff.c Makefile

partdiff-openmp.o: Makefile
	$(CC) -c $(CFLAGS) -fopenmp -o partdiff-openmp.o partdiff.c

partdiff-openmp-zeile.o: Makefile
	$(CC) -c $(CFLAGS) -D ZEILE -fopenmp -o partdiff-openmp-zeile.o partdiff.c

partdiff-openmp-spalte.o: Makefile
	$(CC) -c $(CFLAGS) -D SPALTE -fopenmp -o partdiff-openmp-spalte.o partdiff.c

partdiff-openmp-element.o: Makefile
	$(CC) -c $(CFLAGS) -D ELEMENT -fopenmp -o partdiff-openmp-element.o partdiff.c

askparams.o: askparams.c Makefile

clean:
	$(RM) *.o *~
	$(RM) $(TGTS)

