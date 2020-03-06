#include<uGnix.h>
#include<limits.h>
#include<assert.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

gsl_rng * r;

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


static chromosome* getChrPtr(int chr, chrsample* chrom)
{
  chromosome* tmp = chrom->chrHead;
  int i = 0;
  while( i <= chr )
    {
      assert(tmp != NULL);
      if(i==chr)
	break;
      tmp = tmp->next;
      i++;
    }
  return(tmp);
}

static unsigned int unionAnc(unsigned int anc1, unsigned int anc2)
{
  return(anc1 | anc2);
}

static void delete_anc(ancestry* head)
{
  struct ancestry* tmp;
   while (head != NULL)
    {
       tmp = head;
       head = head->next;
       free(tmp);
    }
}

static void delete_chrom(chromosome* chrptr, chrsample* chrom)
{
  chromosome* tempChrom;
  chromosome* currChrom = chrom->chrHead;
  
  if(currChrom == chrptr)
    {
      tempChrom = currChrom;
      currChrom = currChrom->next;
      chrom->chrHead = currChrom;
    }
  else
    {
      while(chrptr != currChrom->next)
	currChrom = currChrom->next;
      tempChrom = currChrom->next;
      currChrom->next = currChrom->next->next;
    }
  delete_anc(tempChrom->anc);
  free(tempChrom);
}

static void delete_sample(chromosome* head)
{
  chromosome* tmp;
  while (head != NULL)
    {
      tmp = head;
      head = head->next;
      delete_anc(tmp->anc);
      free(tmp);
    }
}

static double totalAncLength(const chrsample* chrom)
{
  chromosome* currChrom = chrom->chrHead;
  ancestry* tmp_anc; 
  double lastPosition = 0;
  double totLength = 0;
  while(currChrom != NULL)
    {
      tmp_anc = currChrom->anc;
      lastPosition=0;
      while(tmp_anc != NULL)
	{
	  if(tmp_anc->abits)
	      totLength += tmp_anc->position - lastPosition;
	  lastPosition = tmp_anc->position;
	  tmp_anc = tmp_anc->next;
	}
      currChrom = currChrom->next;
    }
  return(totLength);
}

void getRecEvent(chrsample* chrom, double eventPos, recombination_event* recEv)
{
  recEv->location = 0.8;
  recEv->chrom = getChrPtr(4, chrom);
}


void recombination(unsigned int* noChrom, recombination_event recEv, chrsample* chrom)
{
  chromosome* chrtmp;
  chrtmp = recEv.chrom;
  ancestry* tmp = chrtmp->anc;
  ancestry* currAnc = NULL;
  chromosome* newLeft = malloc(sizeof(chromosome));
  chromosome* newRight = malloc(sizeof(chromosome));
  newLeft->next = NULL;
  newRight->next = NULL;

  newRight->anc = malloc(sizeof(ancestry));
  newRight->anc->abits = 0;
  newRight->anc->position = recEv.location;
  newRight->anc->next = NULL;
  currAnc = newRight->anc;
  while( tmp != NULL )
    {
      if(recEv.location < tmp->position)
	{
	  currAnc->next = malloc(sizeof(ancestry));
	  currAnc = currAnc->next;
	  currAnc->next = NULL;
	  currAnc->abits = tmp->abits;
	  currAnc->position = tmp->position;
	}
      tmp = tmp->next;
    }

  newLeft->anc = malloc(sizeof(ancestry));
  newLeft->anc->next = NULL;
  tmp = chrtmp->anc;
  currAnc = newLeft->anc;
  int atHead = 1;
  while( tmp->position < recEv.location )
    {
      if( !atHead )
	{
	  currAnc->next = malloc(sizeof(ancestry));
	  currAnc = currAnc->next;
	  currAnc->next = NULL;
	}
      currAnc->abits = tmp->abits;
      currAnc->position = tmp->position;
      tmp = tmp->next;
      atHead = 0;
    }
  if( !atHead )
    {
      currAnc->next = malloc(sizeof(ancestry));
      currAnc = currAnc->next;
      currAnc->next = NULL;
    }
  currAnc->abits = tmp->abits;
  currAnc->position = recEv.location;
  currAnc->next = malloc(sizeof(ancestry));
  currAnc = currAnc->next;
  currAnc->abits = 0;
  currAnc->position = 1.0;
  currAnc->next = NULL;
    
  delete_chrom(chrtmp, chrom);
  chrtmp = chrom->chrHead;
  while(chrtmp->next != NULL)
    chrtmp = chrtmp->next;
  
  currAnc = newLeft->anc;
  unsigned int sumAnc=0;
  while( currAnc != NULL)
    {
      sumAnc += currAnc->abits;
      currAnc = currAnc->next;
    }
  if( sumAnc != 0 ) // check that chromosome is ancestral to sample -- otherwise discard
    {
      chrtmp->next = newLeft;
      chrtmp = chrtmp->next;
      *noChrom = *noChrom + 1;
    }
  else
    {
      delete_anc(newLeft->anc);
      free(newLeft);
    }

  currAnc = newRight->anc;
  sumAnc=0;
  while( currAnc != NULL)
    {
      sumAnc += currAnc->abits;
      currAnc = currAnc->next;
    }
  if( sumAnc != 0 ) // check that chromosome is ancestral to sample -- otherwise discard
    {
      chrtmp->next = newRight;
      *noChrom = *noChrom + 1;
    }
    else
    {
      delete_anc(newRight->anc);
      free(newRight);
    }
  *noChrom = *noChrom - 1;
}

chromosome* mergeChr(chromosome* ptrchr1, chromosome* ptrchr2)
{
  double epsilon = 0.000001;
  chromosome* commonAnc = malloc(sizeof(chromosome));
  commonAnc->next = NULL;
  commonAnc->anc = malloc(sizeof(ancestry));
  commonAnc->anc->abits=0;
  commonAnc->anc->position=0;
  commonAnc->anc->next = NULL; 
  ancestry* tmp = commonAnc->anc;
  ancestry* anc1 = ptrchr1->anc;
  ancestry* anc2 = ptrchr2->anc;
  
  while((anc1 != NULL)&&(anc2 != NULL))
    {
      tmp->abits = unionAnc(anc1->abits,anc2->abits);
      if((anc1->position - anc2->position) > epsilon)
	{
	  tmp->position = anc2->position;
	  anc2 = anc2->next;
	}
      else
	{
	  if((anc2->position - anc1->position) > epsilon)
	    {
	      tmp->position = anc1->position;
	      anc1 = anc1->next;
	    }
	  else
	    {
	      tmp->position = anc1->position;
	      anc2 = anc2->next;
	      anc1 = anc1->next;
	    }
	}
       if((anc1 != NULL)&&(anc2 != NULL))
	 {
	   tmp->next = malloc(sizeof(ancestry));
	   tmp->next->next = NULL;
	   tmp->next->abits=0;
	   tmp->next->position=0;
	   tmp = tmp->next;
	 }
    } 
  return(commonAnc);
} 

static void combineIdentAdjAncSegs(chromosome *ptrchr)
{
  ancestry* tmp;
  ancestry* tmp_del;
  tmp = ptrchr->anc;
  while(tmp->next != NULL)
    {
      if(tmp->abits == tmp->next->abits)
	{
	  tmp_del = tmp->next;
	  tmp->position = tmp->next->position;
	  tmp->next = tmp->next->next;
	  free(tmp_del);
	}
      if(tmp->next != NULL)
	tmp = tmp->next;
    }
}

static void coalescence(unsigned int* noChrom, chrsample* chrom)
{

  int chr1; 
  int chr2;
  if(*noChrom > 2)
    {
      chr1 = gsl_rng_uniform_int(r, *noChrom - 1);
      chr2 = gsl_rng_uniform_int(r, *noChrom - 2);
      if(chr2 >= chr1)
	chr2++;
    }
  else
    {
      chr1 = 0;
      chr2 = 1;
    }
  *noChrom = *noChrom - 1;
  chromosome* tmp;
  chromosome* ptrchr1 = NULL;
  chromosome* ptrchr2 = NULL;
  chromosome* commonAnc = NULL;
  ptrchr1 = getChrPtr(chr1, chrom);
  ptrchr2 = getChrPtr(chr2, chrom);
  commonAnc = mergeChr(ptrchr1, ptrchr2);
  combineIdentAdjAncSegs(commonAnc);
  delete_chrom(ptrchr1,chrom);
  delete_chrom(ptrchr2,chrom);
  if(*noChrom > 1)
    {
      tmp = chrom->chrHead;
      while(tmp->next != NULL)
	tmp = tmp->next;
      tmp->next = commonAnc;
    }
  else
    chrom->chrHead = commonAnc;
}

static chrsample* create_sample(int noChrom)
{
  chrsample* chromSample = malloc(sizeof(chrsample));
  chromSample->chrHead = NULL;
  chromosome* currentChrom;
  chromosome* newChrom;

  // create linked list of n sampled chromosomes
  for(int i=0; i<noChrom; i++)
    {
      if(chromSample->chrHead == NULL)
	{
	  chromSample->chrHead = malloc(sizeof(chromosome));
	  chromSample->chrHead->next = NULL;
	  chromSample->chrHead->anc = malloc(sizeof(ancestry));
	  chromSample->chrHead->anc->next = NULL;
	  chromSample->chrHead->anc->abits = 1;
	  chromSample->chrHead->anc->position = 1.0;
	  currentChrom = chromSample->chrHead;
	}
      else
	{
	  newChrom = malloc(sizeof(chromosome));
	  newChrom->next = NULL;
	  newChrom->anc = malloc(sizeof(ancestry));
	  newChrom->anc->next = NULL;
	  newChrom->anc->abits = 1;
	  newChrom->anc->abits <<= i;
	  newChrom->anc->position = 1.0;
	  currentChrom->next = newChrom;
	  currentChrom = newChrom;
	}
    }
  return(chromSample);
}

char version[] = "coalsim";

int main()
{
  const gsl_rng_type * T;

 /* create a generator chosen by the
    environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();

  // gsl_rng_default_seed = 45567; 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  recombination_event recombEvent;
  unsigned int noChrom=20;
  unsigned int noSamples=20;
  //  unsigned int MRCA=0;
  chromosome* currentChrom = NULL;
  chrsample* chromSample = create_sample(noChrom);
  ancestry* tmp = NULL;
  int currChr=0;
  currentChrom = chromSample->chrHead;
  //  recombination(3, 0.6, chromSample);

  // coalescence(&noChrom, 19, 20, chromSample);
  // 
  // coalescence(&noChrom, 19, 20, chromSample);
  // coalescence(&noChrom, 18, 19, chromSample);
  int noRec=0;
      /*  while(noChrom > 1)
  {
    if(gsl_rng_uniform_pos(r) > 0.85)
      coalescence(&noChrom, chromSample);
    else
      {
	recombination(&noChrom, gsl_rng_uniform_pos(r), chromSample);
	noRec++;
      }
      } */

  getRecEvent(chromSample, 0.6, &recombEvent);
  recombination(&noChrom,recombEvent,chromSample);

  // coalescence(&noChrom, chromSample);
  currChr=0;
  currentChrom = chromSample->chrHead; 
  while(currentChrom != NULL)
    {
      printf("\nChr: %d Anc: ",currChr);
      tmp = currentChrom->anc;
      while(tmp != NULL)
	{
	  displayBits(tmp->abits,noSamples);
	  printf(" %lf ",tmp->position);  
	  tmp = tmp->next;
	  } 
      currentChrom = currentChrom->next;
      currChr++;
      }
  printf("\nTotalanclength: %lf",totalAncLength(chromSample));
  printf("\nNo Rec: %d",noRec);
  printf("\nNoChrom: %d\n",noChrom);
  //  printf("U(0,1): %lf\n",gsl_rng_uniform_pos(r));
  // printf("U(0,1): %lf",gsl_ran_flat(r,0,1.0));
  delete_sample(chromSample->chrHead);
  free(chromSample); 
}

