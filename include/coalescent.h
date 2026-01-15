#include<limits.h>
#include<assert.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include "bitarray.h"

#define POS2BASE(X,Y) (long)(ceil((X)*(Y))-1) /* convert POS X in (0,1) to BASE position
				   in sequence of length Y */

/* Global sample size for bitarray allocation */
extern int g_noSamples;

struct tree {
  struct tree* left;
  struct tree* right;
  bitarray* abits;
  double time;
};

struct mutation
{
  double location;
  bitarray* abits;
  double age;
  struct mutation* next;
};
typedef struct mutation mutation;

struct ancestry
{
  double position;
  bitarray* abits;
  struct ancestry* next;
  int is_mrca;  /* cached: 1 if segment has reached MRCA (permanent) */
};
typedef struct ancestry ancestry;

struct chromosome
{
  ancestry* anc;
  double ancLen;  /* cached ancestral length for O(log n) event lookup */
};
typedef struct chromosome chromosome;

typedef struct
{
  chromosome** chrs;    /* array of chromosome pointers for O(1) access */
  int count;            /* current number of chromosomes */
  int capacity;         /* allocated capacity */
  double ancLength;     /* cached total ancestral length */
  double activeAncLength; /* cached active length (excluding MRCA segments) */
  int activeAncValid;   /* 1 if activeAncLength is valid, 0 if needs recalc */
} chrsample;

typedef struct
{
  double location;
  chromosome* chrom;
  int chromIdx;  /* index for O(1) deletion */
} recombination_event;

typedef struct
{
  int chr1;
  int chr2;
} coalescent_pair;

struct coalescent_events {
  chromosome* chr;
  double time;
  struct coalescent_events* next;
};

struct geneTree {
  bitarray* abits;
  double time;
  struct geneTree* next;
};

struct mrca_list {
  double lower_end;
  double upper_end;
  double age;
  struct mrca_list* next;
};

struct mrca_summary {
  double age;
  double length;
  int numInts;
  struct mrca_summary* next;
};

/* O(1) chromosome access via array */
chromosome* getChrPtr(int chr, chrsample* chrom);

/* Union two bitarrays (creates new bitarray with combined bits) */
bitarray* unionAnc(const bitarray *anc1, const bitarray *anc2);

chromosome* copy_chrom(chromosome* sourceChr);

void delete_anc(ancestry* head);

/* O(1) swap-and-pop deletion */
void delete_chrom_idx(int idx, chrsample* chrom);

void delete_sample(chrsample* chrom);

/* Calculate total ancestral length (updates cache) */
double calcAncLength(chrsample* chrom);

/* Calculate ancestral length for a single chromosome */
double calcChromAncLength(const chromosome* chr);

/* Calculate active ancestral length (excluding MRCA segments) for a chromosome */
double calcChromAncLengthActive(const chromosome* chr, const bitarray* mrca);

/* Calculate total ACTIVE ancestral length (excluding segments at MRCA) */
double calcActiveAncLength(chrsample* chrom, const bitarray* mrca);

/* Count active segments and their total length (for progress estimation) */
int countActiveSegments(chrsample* chrom, const bitarray* mrca, double* totalLen);

/* Get cached ancestral length */
double getAncLength(const chrsample* chrom);

/* Update ancLength after coalescence (pass the overlap that was removed) */
void updateAncLengthCoal(chrsample* chrom, double removed);

/* Get active ancestral length with caching.
 * Returns cached value if valid, otherwise recalculates.
 * Cache is invalidated after coalescence (when MRCA segments can change). */
double getActiveAncLength(chrsample* chrom, const bitarray* mrca);

/* Invalidate active ancestral length cache (call after coalescence) */
void invalidateActiveAncLength(chrsample* chrom);

/* Append chromosome to sample (grows array if needed) */
void appendChrom(chrsample* chrom, chromosome* chr);

void getRecEvent(chrsample* chrom, double eventPos, recombination_event* recEv);

/* Get recombination event, skipping MRCA segments */
void getRecEventActive(chrsample* chrom, double eventPos, recombination_event* recEv,
                       const bitarray* mrca);

void recombination(unsigned int* noChrom, recombination_event recEv, chrsample* chrom);

chromosome* mergeChr(chromosome* ptrchr1, chromosome* ptrchr2);

void combineIdentAdjAncSegs(chromosome *ptrchr);

void coalescence(coalescent_pair pair, unsigned int* noChrom, chrsample* chrom);

void updateCoalescentEvents(struct coalescent_events** coalescent_list,
			    struct coalescent_events** coalescent_list_tail,
			    chrsample* chromSample, double totalTime);

struct geneTree* getGeneTree(double lower, double upper, struct coalescent_events* coalescent_list, const bitarray* mrca);

unsigned long long int ipow( unsigned long long int base, int exp);

int TestMRCAForAll(chrsample* chrom, const bitarray* mrca);

chrsample* create_sample(int noChrom);

void getCoalPair(gsl_rng * r, unsigned int noChrom, coalescent_pair* pair);

void addMRCAInterval(struct mrca_list** head, double newlower,
		     double newupper, double newage);

void getMRCAs(struct mrca_list** head, chrsample* chromSample, double totalTime, const bitarray* mrca);

void MRCAStats(struct mrca_list* head, struct mrca_summary* mrca_head, double smalldiff, long chromTotBases,
		  int seqUnits, char* baseUnit, int prn_mrca, int prn_regions);

void getMutEvent(chrsample* chrom, double eventPos, mutation* mutEv, double time);

/* Get mutation event, skipping MRCA segments */
void getMutEventActive(chrsample* chrom, double eventPos, mutation* mutEv, double time,
                       const bitarray* mrca);

void printMutations(mutation* mutation_list, long chromTotBases, int seqUnits,
		    char* baseUnit, int noSamples, const bitarray* mrca);

void printChromosomes(chrsample* chromSample, int noSamples);

long convertToBases(long totBases, int seqUnit, double value);

/* Substitution model types */
typedef enum {
  SUBST_JC69,   /* Jukes-Cantor: equal rates, equal frequencies */
  SUBST_HKY     /* HKY85: transition/transversion bias, variable frequencies */
} subst_model_t;

/* HKY model parameters */
typedef struct {
  double kappa;      /* transition/transversion ratio */
  double pi[4];      /* base frequencies: A, C, G, T */
  double ti_prob[4]; /* precomputed transition probs (normalized) */
  double tv_prob[4][2]; /* precomputed transversion probs per base */
} hky_params_t;

/* Initialize HKY parameters (call once before simulation) */
void initHKYParams(hky_params_t* params, double kappa,
                   double piA, double piC, double piG, double piT);

/* Substitution functions */
char JC69RBase(gsl_rng * r, char currBase);
char HKYRBase(gsl_rng * r, char currBase, const hky_params_t* params);

/* Sequence simulation with model selection */
char** simulateSequences(mutation* mutation_list, int totBases, int noSamples,
                         gsl_rng * r, subst_model_t model, const hky_params_t* hky);

/* VCF output - writes only variable sites (memory efficient) */
void writeVCF(mutation* mutation_list, long chromTotBases, int noSamples,
              const bitarray* mrca, FILE* out, gsl_rng* r);

/*
 * Target regions for sparse coalescent simulation.
 * When target regions are specified, the simulation terminates once
 * all target regions have reached MRCA, rather than waiting for the
 * entire chromosome. This can dramatically speed up simulations for
 * long chromosomes when only specific regions are of interest.
 */
typedef struct {
    double start;    /* Start position (0-1 scale) */
    double end;      /* End position (0-1 scale) */
} target_region;

typedef struct {
    target_region* regions;
    int n_regions;
    int active;      /* 1 if sparse mode, 0 for full simulation */
    double min_start; /* Cached: minimum start across all regions */
    double max_end;   /* Cached: maximum end across all regions */
} target_region_set;

/* Create target region set from parameters:
 * n_regions: number of target regions
 * region_length: length of each region (in 0-1 scale, e.g., 0.0001 for 0.01 cM)
 * first_start: start position of first region
 * spacing: distance between region starts (0 = evenly distributed)
 */
target_region_set* create_target_regions(int n_regions, double region_length,
                                         double first_start, double spacing);

/* Create target region set from file (one region per line: start end) */
target_region_set* load_target_regions(const char* filename);

/* Free target region set */
void free_target_regions(target_region_set* regions);

/* Check if a position (0-1 scale) falls within any target region.
 * Returns 1 if in a target region, 0 otherwise.
 * If regions is NULL or not active, returns 1 (accept all positions).
 */
int position_in_target_regions(double pos, const target_region_set* regions);

/* Test if all target regions have reached MRCA.
 * Returns 1 if all target regions have MRCA, 0 otherwise.
 * If regions is NULL or not active, falls back to TestMRCAForAll behavior.
 */
int TestMRCAForTargetRegions(chrsample* chrom, const bitarray* mrca,
                             const target_region_set* regions);

/* Enable sparse ancestry tracking for target regions.
 * When set, bitarray operations are skipped for segments outside target regions.
 * This dramatically reduces computation for sparse simulations. */
void set_sparse_target_regions(const target_region_set* trs);

/* Check if chromosome has any segment overlapping target regions.
 * Used to prune chromosomes that can't contribute to target genealogy. */
int chromosome_overlaps_targets(const chromosome* chr);

/* Deprecated: use bitarray functions instead */
void getBits(unsigned int value, unsigned int noSamples, unsigned int* result);

/* Use bitarray_is_singleton instead for bitarray types */
int isSingleton(unsigned int x);

void addNode(bitarray* val, double time, struct tree* lroot);

void splitNode(struct tree* lroot);

/* Convert bitarray to chromosome label (returns first set bit index + 1) */
int binaryToChrLabel(const bitarray* ba);

void fillTips(struct tree* lroot);

void printTree(struct tree* lroot, int noSamples, int toScreen, FILE* tree_file);

void printTreeNewick(struct tree* root, int noSamples, int toScreen, FILE* tree_file);

void freeTree(struct tree* node);


