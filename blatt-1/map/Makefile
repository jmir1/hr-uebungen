all: heatmap
	
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb

heatmap: heatmap.o
	gcc $(CFLAGS) -o heatmap heatmap.o
	
heatmap.o: heatmap.c
	gcc $(CFLAGS) -c heatmap.c

clean:
	rm -f heatmap
	rm -f *.o
	rm -f *~
