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

partdiff-openmp: partdiff-openmp.o askparams-openmp.o Makefile
	gcc $(LFLAGS) -fopenmp -DOPENMP -o $@ partdiff-openmp.o askparams-openmp.o -lm

partdiff-openmp-element: partdiff-openmp-element.o askparams-openmp.o Makefile
	gcc $(LFLAGS) -fopenmp -DOPENMP -DELEMENT -o $@ partdiff-openmp-element.o askparams-openmp.o -lm

partdiff-openmp-spalte: partdiff-openmp-spalte.o askparams-openmp.o Makefile
	gcc $(LFLAGS) -fopenmp -DOPENMP -DSPALTE -o $@ partdiff-openmp-spalte.o askparams-openmp.o -lm

partdiff-openmp-zeile: partdiff-openmp-zeile.o askparams-openmp.o Makefile
	gcc $(LFLAGS) -fopenmp -DOPENMP -DZEILE -o $@ partdiff-openmp-zeile.o askparams-openmp.o -lm



partdiff.o: partdiff.c Makefile

partdiff-openmp.o: Makefile
	$(CC) -c $(CFLAGS) -DOPENMP -fopenmp -o partdiff-openmp.o partdiff.c

partdiff-openmp-zeile.o: Makefile
	$(CC) -c $(CFLAGS) -DOPENMP -DZEILE -fopenmp -o partdiff-openmp-zeile.o partdiff.c

partdiff-openmp-spalte.o: Makefile
	$(CC) -c $(CFLAGS) -DOPENMP -DSPALTE -fopenmp -o partdiff-openmp-spalte.o partdiff.c

partdiff-openmp-element.o: Makefile
	$(CC) -c $(CFLAGS) -DOPENMP -DELEMENT -fopenmp -o partdiff-openmp-element.o partdiff.c

askparams.o: askparams.c Makefile

askparams-openmp.o: askparams.c Makefile
	$(CC) -c $(CFLAGS) -DOPENMP -o askparams-openmp.o askparams.c

clean:
	$(RM) *.o *~
	$(RM) $(TGTS)

