#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "vcfassemble.h"

static void print_usage(const char* progname)
{
    fprintf(stderr, "Usage: %s -v <founder.vcf> -p <pedtrans.txt> -o <output.vcf> [-s]\n", progname);
    fprintf(stderr, "\nAssemble sample VCF from founder VCF and pedtrans segments\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -v FILE   Founder VCF file from coalsim\n");
    fprintf(stderr, "  -p FILE   Pedtrans segment output file\n");
    fprintf(stderr, "  -o FILE   Output VCF file\n");
    fprintf(stderr, "  -s        Output samples only (exclude founders)\n");
    fprintf(stderr, "  -h        Show this help message\n");
}

int main(int argc, char** argv)
{
    char* vcf_file = NULL;
    char* pedtrans_file = NULL;
    char* output_file = NULL;
    int samples_only = 0;

    int c;
    while ((c = getopt(argc, argv, "v:p:o:sh")) != -1) {
        switch (c) {
            case 'v':
                vcf_file = optarg;
                break;
            case 'p':
                pedtrans_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 's':
                samples_only = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (!vcf_file || !pedtrans_file || !output_file) {
        fprintf(stderr, "Error: missing required arguments\n");
        print_usage(argv[0]);
        return 1;
    }

    /* Read founder VCF */
    fprintf(stderr, "Reading founder VCF: %s\n", vcf_file);
    founder_vcf* fvcf = read_founder_vcf(vcf_file);
    if (!fvcf) {
        return 1;
    }
    fprintf(stderr, "  %d mutations, %d founder haplotypes, %d chromosomes\n",
            fvcf->n_mutations, fvcf->n_founder_samples, fvcf->n_chromosomes);

    /* Read pedtrans segments */
    fprintf(stderr, "Reading pedtrans segments: %s\n", pedtrans_file);
    pedtrans_segments* pts = read_pedtrans_segments(pedtrans_file);
    if (!pts) {
        free_founder_vcf(fvcf);
        return 1;
    }
    fprintf(stderr, "  %d sample chromosomes, %d unique samples, %d chromosomes\n",
            pts->n_chroms, pts->n_samples, pts->n_chromosomes);

    /* Open output file */
    FILE* out = fopen(output_file, "w");
    if (!out) {
        fprintf(stderr, "Error: cannot open output file %s\n", output_file);
        free_founder_vcf(fvcf);
        free_pedtrans_segments(pts);
        return 1;
    }

    /* Assemble sample VCF */
    fprintf(stderr, "Assembling sample VCF...\n");
    int ret = assemble_sample_vcf(fvcf, pts, out, samples_only);

    fclose(out);
    fprintf(stderr, "Output written to: %s\n", output_file);

    /* Cleanup */
    free_founder_vcf(fvcf);
    free_pedtrans_segments(pts);

    return ret;
}
