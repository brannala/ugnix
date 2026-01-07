#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "pedsim_vcf.h"

static void print_usage(const char* progname)
{
    fprintf(stderr, "Usage: %s [options] -o <output.vcf>\n", progname);
    fprintf(stderr, "\nIntegrated pedigree-coalescent simulation with VCF output\n\n");
    fprintf(stderr, "Pedigree options:\n");
    fprintf(stderr, "  -g NUM    Number of generations to trace back (default: 5)\n");
    fprintf(stderr, "  -n NUM    Sample size - individuals at present (default: 10)\n");
    fprintf(stderr, "  -P NUM    Pedigree population size (default: 1000)\n");
    fprintf(stderr, "\nCoalescent options:\n");
    fprintf(stderr, "  -N NUM    Effective population size for coalescent (default: 10000)\n");
    fprintf(stderr, "  -r NUM    Recombination rate in cM (default: 1.0)\n");
    fprintf(stderr, "  -m NUM    Mutation rate (default: 0.5)\n");
    fprintf(stderr, "  -s NUM    Random seed (default: 0 = random)\n");
    fprintf(stderr, "\nOutput options:\n");
    fprintf(stderr, "  -o FILE   Output VCF file (required)\n");
    fprintf(stderr, "  -S        Output samples only (exclude founders)\n");
    fprintf(stderr, "  -k        Keep temporary files\n");
    fprintf(stderr, "  -v        Verbose output\n");
    fprintf(stderr, "  -h        Show this help message\n");
    fprintf(stderr, "\nExample:\n");
    fprintf(stderr, "  %s -g 5 -n 10 -N 10000 -r 1.0 -m 0.5 -o output.vcf\n", progname);
}

int main(int argc, char** argv)
{
    pedsim_vcf_params params = {
        .n_generations = 5,
        .sample_size = 10,
        .ped_pop_size = 1000,
        .coal_pop_size = 10000,
        .rec_rate = 1.0,
        .mut_rate = 0.5,
        .seed = 0,
        .output_vcf = NULL,
        .samples_only = 0,
        .keep_temp = 0,
        .verbose = 0
    };

    int c;
    while ((c = getopt(argc, argv, "g:n:P:N:r:m:s:o:Skvh")) != -1) {
        switch (c) {
            case 'g':
                params.n_generations = atoi(optarg);
                break;
            case 'n':
                params.sample_size = atoi(optarg);
                break;
            case 'P':
                params.ped_pop_size = atof(optarg);
                break;
            case 'N':
                params.coal_pop_size = atof(optarg);
                break;
            case 'r':
                params.rec_rate = atof(optarg);
                break;
            case 'm':
                params.mut_rate = atof(optarg);
                break;
            case 's':
                params.seed = atoi(optarg);
                break;
            case 'o':
                params.output_vcf = optarg;
                break;
            case 'S':
                params.samples_only = 1;
                break;
            case 'k':
                params.keep_temp = 1;
                break;
            case 'v':
                params.verbose = 1;
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (!params.output_vcf) {
        fprintf(stderr, "Error: output VCF file required (-o)\n");
        print_usage(argv[0]);
        return 1;
    }

    fprintf(stderr, "Pedigree-Coalescent VCF Simulation\n");
    fprintf(stderr, "===================================\n");
    fprintf(stderr, "Pedigree: %d generations, %d samples, N_ped=%.0f\n",
            params.n_generations, params.sample_size, params.ped_pop_size);
    fprintf(stderr, "Coalescent: N=%.0f, r=%.2f cM, m=%.2f\n",
            params.coal_pop_size, params.rec_rate, params.mut_rate);

    return run_pedsim_vcf_pipeline(&params);
}
