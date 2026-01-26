LDFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I include
TESTFLAGS = -I testing
CC = gcc
vpath %.c src testing
vpath %.h include testing
CFLAGS = -Wall -O3
PROFILE = -g

# build programs

all: hwe-dis kinship gsum het coalsim coalsim_msc pedtrans pedsim pedsim_multipop seqassemble pedsim_seq vcfassemble pedsim_vcf pedsim_vcf_multipop sample
hwe-dis: hwe-dis.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) hwe-dis.o uGnix.o -lglib-2.0 -lm -o hwe-dis
kinship: kinship.o data.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) kinship.o data.o uGnix.o -lglib-2.0 -lm -o kinship
het: het.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) het.o uGnix.o -lglib-2.0 -lm -o het
gsum: gsum.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) gsum.o uGnix.o -lglib-2.0 -lm -o gsum
coalsim: coalsim.o coalescent.o fenwick.o bitarray.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) coalsim.o coalescent.o fenwick.o bitarray.o uGnix.o -lglib-2.0 -lm -lgsl -lgslcblas -o coalsim
coalsim_msc: msc_main.o msc.o species_tree.o coalescent.o fenwick.o bitarray.o uGnix.o -lglib-2.0
	$(CC) $(PROFILE) msc_main.o msc.o species_tree.o coalescent.o fenwick.o bitarray.o uGnix.o -lglib-2.0 -lm -lgsl -lgslcblas -o coalsim_msc
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
coalsim.o: coalsim.c uGnix.h coalescent.h bitarray.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
msc_main.o: msc_main.c msc.h species_tree.h coalescent.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
msc.o: msc.c msc.h species_tree.h coalescent.h bitarray.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
species_tree.o: species_tree.c species_tree.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
coalescent.o: coalescent.c uGnix.h coalescent.h bitarray.h fenwick.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
fenwick.o: fenwick.c fenwick.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
bitarray.o: bitarray.c bitarray.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedtrans: pedtrans_main.o pedtrans.o -lglib-2.0 -lgsl -lgslcblas
	$(CC) $(PROFILE) pedtrans_main.o pedtrans.o -lglib-2.0 -lm -lgsl -lgslcblas -o pedtrans
pedtrans_main.o: pedtrans_main.c pedtrans.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedtrans.o: pedtrans.c pedtrans.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim: pedsim_main.o pedsim.o -lgsl -lgslcblas
	$(CC) $(PROFILE) pedsim_main.o pedsim.o -lm -lgsl -lgslcblas -o pedsim
pedsim_main.o: pedsim_main.c pedsim.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim.o: pedsim.c pedsim.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_multipop: pedsim_multipop_main.o pedsim_multipop.o -lgsl -lgslcblas
	$(CC) $(PROFILE) pedsim_multipop_main.o pedsim_multipop.o -lm -lgsl -lgslcblas -o pedsim_multipop
pedsim_multipop_main.o: pedsim_multipop_main.c pedsim_multipop.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_multipop.o: pedsim_multipop.c pedsim_multipop.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
seqassemble: seqassemble_main.o seqassemble.o
	$(CC) $(PROFILE) seqassemble_main.o seqassemble.o -o seqassemble
seqassemble_main.o: seqassemble_main.c seqassemble.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
seqassemble.o: seqassemble.c seqassemble.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_seq: pedsim_seq_main.o pedsim_seq.o
	$(CC) $(PROFILE) pedsim_seq_main.o pedsim_seq.o -o pedsim_seq
pedsim_seq_main.o: pedsim_seq_main.c pedsim_seq.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_seq.o: pedsim_seq.c pedsim_seq.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
vcfassemble: vcfassemble_main.o vcfassemble.o
	$(CC) $(PROFILE) vcfassemble_main.o vcfassemble.o -o vcfassemble
vcfassemble_main.o: vcfassemble_main.c vcfassemble.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
vcfassemble.o: vcfassemble.c vcfassemble.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_vcf: pedsim_vcf_main.o pedsim_vcf.o
	$(CC) $(PROFILE) pedsim_vcf_main.o pedsim_vcf.o -o pedsim_vcf
pedsim_vcf_main.o: pedsim_vcf_main.c pedsim_vcf.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_vcf.o: pedsim_vcf.c pedsim_vcf.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_vcf_multipop: pedsim_vcf_multipop_main.o pedsim_vcf_multipop.o
	$(CC) $(PROFILE) pedsim_vcf_multipop_main.o pedsim_vcf_multipop.o -o pedsim_vcf_multipop
pedsim_vcf_multipop_main.o: pedsim_vcf_multipop_main.c pedsim_vcf_multipop.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
pedsim_vcf_multipop.o: pedsim_vcf_multipop.c pedsim_vcf_multipop.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
sample: sample_main.o sample.o -lgsl -lgslcblas -lglib-2.0
	$(CC) $(PROFILE) sample_main.o sample.o -lm -lgsl -lgslcblas -lglib-2.0 -o sample
sample_main.o: sample_main.c sample.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
sample.o: sample.c sample.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
uGnix.o: uGnix.c
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) -c $<
clean:
	$(RM) gsum het coalsim coalsim_msc test_ugnix test_coalescent test_msc runtests kinship hwe-dis pedtrans test_pedtrans pedsim pedsim_multipop seqassemble pedsim_seq vcfassemble pedsim_vcf pedsim_vcf_multipop sample
	$(RM) gsum.o uGnix.o het.o coalsim.o coalescent.o fenwick.o bitarray.o msc.o msc_main.o species_tree.o test_ugnix.o test_coalescent.o test_msc.o unity.o kinship.o data.o hwe-dis.o pedtrans.o pedtrans_main.o test_pedtrans.o pedsim.o pedsim_main.o pedsim_multipop.o pedsim_multipop_main.o seqassemble.o seqassemble_main.o pedsim_seq.o pedsim_seq_main.o vcfassemble.o vcfassemble_main.o pedsim_vcf.o pedsim_vcf_main.o pedsim_vcf_multipop.o pedsim_vcf_multipop_main.o sample.o sample_main.o
tidy:
	$(RM) *.o

# build test suite

tests: test_ugnix test_coalescent test_pedtrans test_msc runtests
test_ugnix: test_ugnix.o uGnix.o unity.o -lglib-2.0
	$(CC) $(PROFILE) test_ugnix.o uGnix.o unity.o -lglib-2.0 -lm -o test_ugnix
test_ugnix.o: test_ugnix.c uGnix.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
unity.o: unity.c unity.h unity_internals.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
test_coalescent: test_coalescent.o coalescent.o fenwick.o bitarray.o uGnix.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas
	$(CC) $(PROFILE) test_coalescent.o coalescent.o fenwick.o bitarray.o uGnix.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas -o test_coalescent
test_coalescent.o: test_coalescent.c coalescent.h bitarray.h fenwick.h uGnix.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
test_pedtrans: test_pedtrans.o pedtrans.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas
	$(CC) $(PROFILE) test_pedtrans.o pedtrans.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas -o test_pedtrans
test_pedtrans.o: test_pedtrans.c pedtrans.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
test_msc: test_msc.o msc.o species_tree.o coalescent.o fenwick.o bitarray.o uGnix.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas
	$(CC) $(PROFILE) test_msc.o msc.o species_tree.o coalescent.o fenwick.o bitarray.o uGnix.o unity.o -lglib-2.0 -lm -lgsl -lgslcblas -o test_msc
test_msc.o: test_msc.c msc.h species_tree.h coalescent.h bitarray.h fenwick.h uGnix.h unity.h
	$(CC) $(PROFILE) $(CFLAGS) $(LDFLAGS) $(TESTFLAGS) -c $<
runtests:
	@eval $$(echo "#!/bin/bash" > runtests)
	@eval $$(echo "echo \"\nRunning tests on ugnix.c ...\"" >> runtests)
	@eval $$(echo "./test_ugnix" >> runtests)
	@eval $$(echo "echo \"\nRunning tests on coalescent.c ...\"" >> runtests)
	@eval $$(echo "./test_coalescent" >>runtests)
	@eval $$(echo "echo \"\nRunning tests on pedtrans.c ...\"" >> runtests)
	@eval $$(echo "./test_pedtrans" >>runtests)
	@eval $$(echo "echo \"\nRunning tests on msc.c ...\"" >> runtests)
	@eval $$(echo "./test_msc" >>runtests)
	@eval $$(chmod +x runtests)
testsclean:
	$(RM) test_ugnix
	$(RM) test_coalescent
	$(RM) test_pedtrans
	$(RM) test_msc
	$(RM) runtests
	$(RM) test_ugnix.o unity.o uGnix.o test_coalescent.o coalescent.o test_pedtrans.o pedtrans.o test_msc.o msc.o species_tree.o

