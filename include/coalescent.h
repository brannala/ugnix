#include<limits.h>
#include<assert.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

struct ancestry
{
  double position;
  unsigned int abits;
  struct ancestry* next;
};
typedef struct ancestry ancestry;

struct chromosome
{
  ancestry* anc;
  struct chromosome* next;
};
typedef struct chromosome chromosome;

typedef struct
{
  chromosome* chrHead;
} chrsample;

typedef struct
{
  double location;
  chromosome* chrom;
} recombination_event;

typedef struct
{
  int chr1;
  int chr2;
} coalescent_pair;

chromosome* getChrPtr(int chr, chrsample* chrom);

unsigned int unionAnc(unsigned int anc1, unsigned int anc2);

void delete_anc(ancestry* head);

void delete_chrom(chromosome* chrptr, chrsample* chrom);

void delete_sample(chromosome* head);

double totalAncLength(const chrsample* chrom);

void getRecEvent(chrsample* chrom, double eventPos, recombination_event* recEv);

void recombination(unsigned int* noChrom, recombination_event recEv, chrsample* chrom);

chromosome* mergeChr(chromosome* ptrchr1, chromosome* ptrchr2);

void combineIdentAdjAncSegs(chromosome *ptrchr);

void coalescence(coalescent_pair pair, unsigned int* noChrom, chrsample* chrom);

unsigned long long int ipow( unsigned long long int base, int exp);

int TestMRCAForAll(chrsample* chrom, unsigned int mrca);

chrsample* create_sample(int noChrom);

void getCoalPair(gsl_rng * r, unsigned int noChrom, coalescent_pair* pair);
