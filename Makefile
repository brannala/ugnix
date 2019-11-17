LDFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include
CC = gcc
CFLAGS = -Wall
gsum: gsum.o uGnix.o -lglib-2.0
	$(CC)  gsum.o uGnix.o -lglib-2.0 -o gsum
gsum.o: gsum.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c gsum.c 
uGnix.o: uGnix.c
	$(CC) $(CFLAGS) $(LDFLAGS) -c uGnix.c
clean:
	$(RM) gsum
	$(RM) gsum.o uGnix.o
