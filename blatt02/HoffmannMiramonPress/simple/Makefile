all: simple
	
CFLAGS = -std=c11 -pedantic -Wall -Wextra -Og -ggdb -gdwarf-4

simple: simple.o
	gcc $(CFLAGS) -o simple simple.o
	
simple.o: simple.c
	gcc $(CFLAGS) -c simple.c

clean:
	rm -f simple
	rm -f *.o
	rm -f *~
