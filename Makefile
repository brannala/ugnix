LDFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I include
TESTFLAGS = -I testing
CC = gcc
vpath %.c src testing
vpath %.h include testing
CFLAGS = -Wall -O3
PROFILE = -g

# build programs

all: hwe-dis kinship gsum het coalsim
hwe-dis: hwe-dis.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) hwe-dis.o uGnix.o -lglib-2.0 -lm -o hwe-dis
kinship: kinship.o data.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) kinship.o data.o uGnix.o -lglib-2.0 -lm -o kinship
het: het.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) het.o uGnix.o -lglib-2.0 -lm -o het
gsum: gsum.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) gsum.o uGnix.o -lglib-2.0 -lm -o gsum
coalsim: coalsim.o coalescent.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) coalsim.o coalescent.o uGnix.o -lglib-2.0 -lm -lgsl -lgslcblas -o coalsim
kinship.o: kinship.c kinship_data.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
data.o: data.c kinship_data.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
hwe-dis.o: hwe-dis.c uGnix.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
het.o: het.c uGnix.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
gsum.o: gsum.c uGnix.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
coalsim.o: coalsim.c uGnix.h coalescent.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
coalescent.o: coalescent.c uGnix.h coalescent.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
uGnix.o: uGnix.c
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
clean:
	$(RM) gsum het coalsim test_ugnix test_coalescent runtests kinship hwe-dis
	$(RM) gsum.o uGnix.o het.o coalsim.o coalescent.o test_ugnix.o test_coalescent.o unity.o kinship.o data.o hwe-dis.o
tidy:
	$(RM) *.o

# build test suite

tests: test_ugnix test_coalescent runtests
test_ugnix: test_ugnix.o uGnix.o unity.o -lglib-2.0
	$(CC) $(PROFILE) test_ugnix.o uGnix.o unity.o -lglib-2.0 -lm -o test_ugnix
test_ugnix.o: test_ugnix.c uGnix.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
unity.o: unity.c unity.h unity_internals.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
test_coalescent: test_coalescent.o coalescent.o uGnix.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas
	$(CC) $(PROFILE) test_coalescent.o coalescent.o uGnix.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas -o test_coalescent
test_coalescent.o: test_coalescent.c coalescent.h uGnix.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
runtests:
	@eval $$(echo "#!/bin/bash" > runtests)
	@eval $$(echo "echo \"\nRunning tests on ugnix.c ...\"" >> runtests)
	@eval $$(echo "./test_ugnix" >> runtests)
	@eval $$(echo "echo \"\nRunning tests on coalescent.c ...\"" >> runtests)
	@eval $$(echo "./test_coalescent" >>runtests)
	@eval $$(chmod +x runtests) 	
testsclean:
	$(RM) test_ugnix
	$(RM) test_coalescent
	$(RM) runtests
	$(RM) test_ugnix.o unity.o uGnix.o test_coalescent.o coalescent.o

