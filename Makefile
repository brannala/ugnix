LDFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include
CC = gcc
CFLAGS = -Wall
PROFILE = -g
all: gsum het
het: het.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) het.o uGnix.o -lglib-2.0 -lm -o het
gsum: gsum.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) gsum.o uGnix.o -lglib-2.0 -lm -o gsum
gsum.o: gsum.c
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c gsum.c
het.o: het.c
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c het.c
uGnix.o: uGnix.c
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c uGnix.c
clean:
	$(RM) gsum het
	$(RM) gsum.o uGnix.o het.o
