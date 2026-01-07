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

/* Get cached ancestral length */
double getAncLength(const chrsample* chrom);

/* Update ancLength after coalescence (pass the overlap that was removed) */
void updateAncLengthCoal(chrsample* chrom, double removed);

/* Append chromosome to sample (grows array if needed) */
void appendChrom(chrsample* chrom, chromosome* chr);

void getRecEvent(chrsample* chrom, double eventPos, recombination_event* recEv);

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


