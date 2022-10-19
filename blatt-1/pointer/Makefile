all: pointer
	
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb

pointer: pointer.o
	gcc $(CFLAGS) -o pointer pointer.o
	
pointer.o: pointer.c
	gcc $(CFLAGS) -c pointer.c

clean:
	rm -f pointer
	rm -f *.o
	rm -f *~
