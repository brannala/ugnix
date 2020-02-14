LDFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I include
TESTFLAGS = -I testing
CC = gcc
vpath %.c src testing
vpath %.h include testing
CFLAGS = -Wall
PROFILE = -g

# build programs

all: gsum het
het: het.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) het.o uGnix.o -lglib-2.0 -lm -o het
gsum: gsum.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) gsum.o uGnix.o -lglib-2.0 -lm -o gsum
coalsim: coalsim.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) coalsim.o uGnix.o -lglib-2.0 -lm -o coalsim
het.o: het.c uGnix.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
gsum.o: gsum.c uGnix.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
coalsim.o: coalsim.c uGnix.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
uGnix.o: uGnix.c
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
clean:
	$(RM) gsum het coalsim test_ugnix
	$(RM) gsum.o uGnix.o het.o coalsim.o test_ugnix.o unity.o
tidy:
	$(RM) *.o

# build test suite

tests: test_ugnix
test_ugnix: test_ugnix.o uGnix.o unity.o -lglib-2.0
	$(CC) $(PROFILE) test_ugnix.o uGnix.o unity.o -lglib-2.0 -lm -o test_ugnix
test_ugnix.o: test_ugnix.c uGnix.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
unity.o: unity.c unity.h unity_internals.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
testsclean:
	$(RM) test_ugnix
	$(RM) test_ugnix.o unity.o uGnix.o
