CC = mpicc
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb
LFLAGS = $(CFLAGS)
TGTS = circle

all: $(TGTS)

$(TGTS): %: %.c
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) $(TGTS)
