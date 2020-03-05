#include<uGnix.h>
#include<limits.h>
#include<assert.h>

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


void recombination(int chr, double recLoc, chrsample* chrom)
{
  chromosome* chrtmp;
  chrtmp = getChrPtr(chr, chrom);
  ancestry* tmp = chrtmp->anc;
  ancestry* currAnc = NULL;
  chromosome* newLeft = malloc(sizeof(chromosome));
  chromosome* newRight = malloc(sizeof(chromosome));
  newLeft->next = NULL;
  newRight->next = NULL;

  newRight->anc = malloc(sizeof(ancestry));
  newRight->anc->abits = 0;
  newRight->anc->position = recLoc;
  newRight->anc->next = NULL;
  currAnc = newRight->anc;
  while( tmp != NULL )
    {
      if(recLoc < tmp->position)
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
  while( tmp->position < recLoc )
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
  currAnc->position = recLoc;
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
    chrtmp->next = newRight;
    else
    {
      delete_anc(newRight->anc);
      free(newRight);
    }

}

chromosome* mergeChr(chromosome* ptrchr1, chromosome* ptrchr2)
{
  double epsilon = 0.01;
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

static void coalescence(unsigned int* noChrom, int chr1, int chr2, chrsample* chrom)
{
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
  tmp = chrom->chrHead;
  while(tmp->next != NULL)
    tmp = tmp->next;
  tmp->next = commonAnc;
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
  unsigned int noChrom=20;
  unsigned int noSamples=20;
  //  unsigned int MRCA=0;
  chromosome* currentChrom = NULL;
  chrsample* chromSample = create_sample(noChrom);
  ancestry* tmp = NULL;
  int currChr=0;
  currentChrom = chromSample->chrHead;
  while(currentChrom != NULL)
    {
      printf("Chr: %d",currChr);
      printf(" Segment boundary: %lf Ancestry Bits: ",currentChrom->anc->position);  
      displayBits(currentChrom->anc->abits,noSamples);
      currentChrom = currentChrom->next;
      currChr++;
      } 
  //  delete_chrom(getChrPtr(0, chromSample), chromSample);
  currChr=0;
  currentChrom = chromSample->chrHead;
  while(currentChrom != NULL)
    {
      printf("Chr: %d",currChr);
      printf(" Segment boundary: %lf Ancestry Bits: ",currentChrom->anc->position);  
      displayBits(currentChrom->anc->abits,noSamples);
      currentChrom = currentChrom->next;
      currChr++;
      } 
  //  recombination(3, 0.6, chromSample);

  // coalescence(&noChrom, 19, 20, chromSample);
  // recombination(10, 0.4, chromSample);
  // coalescence(&noChrom, 19, 20, chromSample);
  // coalescence(&noChrom, 18, 19, chromSample);
  // recombination(17, 0.8, chromSample);
  // coalescence(&noChrom, 19, 20, chromSample);
  coalescence(&noChrom, 1, 2, chromSample);
  coalescence(&noChrom, 17, 18, chromSample);
  coalescence(&noChrom, 0, 1, chromSample);
  coalescence(&noChrom, 15, 16, chromSample);
  recombination(15, 0.8, chromSample);
  coalescence(&noChrom, 16, 1, chromSample);
  coalescence(&noChrom, 14, 15, chromSample);
  recombination(7, 0.8, chromSample);
  coalescence(&noChrom, 3, 15, chromSample);
  coalescence(&noChrom, 13, 14, chromSample);
  currChr=0;
  currentChrom = chromSample->chrHead; 
  while(currentChrom != NULL)
    {
      printf("Chr: %d Ancestry segments: ",currChr);
      tmp = currentChrom->anc;
      while(tmp != NULL)
	{
	  printf(" %lf ",tmp->position);  
	  displayBits(tmp->abits,noSamples);
	  tmp = tmp->next;
	  } 
      currentChrom = currentChrom->next;
      currChr++;
      }
  printf("Totalanclength: %lf",totalAncLength(chromSample));
  delete_sample(chromSample->chrHead);
  free(chromSample); 
}

