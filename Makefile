LDFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include
CC = gcc
CFLAGS = -Wall -g
all: gsum het
het: het.o uGnix.o -lglib-2.0
	$(CC) het.o uGnix.o -lglib-2.0 -lm -o het
gsum: gsum.o uGnix.o -lglib-2.0
	$(CC)  gsum.o uGnix.o -lglib-2.0 -lm -o gsum
gsum.o: gsum.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c gsum.c
het.o: het.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c het.c
uGnix.o: uGnix.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c uGnix.c
clean:
	$(RM) gsum
	$(RM) gsum.o uGnix.o het.o
