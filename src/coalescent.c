#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

chromosome* getChrPtr(int chr, chrsample* chrom)
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

unsigned int unionAnc(unsigned int anc1, unsigned int anc2)
{
  return(anc1 | anc2);
}

void delete_anc(ancestry* head)
{
  struct ancestry* tmp;
   while (head != NULL)
    {
       tmp = head;
       head = head->next;
       free(tmp);
    }
}

void delete_chrom(chromosome* chrptr, chrsample* chrom)
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

void delete_sample(chromosome* head)
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

double totalAncLength(const chrsample* chrom)
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

/* 
eventPos is the absolute position of rec event on cumulative ancestral 
chromosome material (ancLength) -> getRecEvent finds the recombinant chromosome 
and relative position of recombination event on that chromosome modifies recEv 
*/

void getRecEvent(chrsample* chrom, double eventPos, recombination_event* recEv)
{
  chromosome* currChrom = chrom->chrHead;
  ancestry* tmp_anc; 
  double lastPosition = 0;
  double currLength = 0;
  unsigned int foundPosition=0;
  while((currChrom != NULL)&&(!foundPosition))
    {
      tmp_anc = currChrom->anc;
      lastPosition=0;
      double lastLength = currLength;
      double nonAncestralLength=0;
      while((tmp_anc != NULL)&&(!foundPosition))
	{
	  if(tmp_anc->abits)
	      currLength += tmp_anc->position - lastPosition;
	  else
	    nonAncestralLength += tmp_anc->position - lastPosition;
	  if(currLength > eventPos)
	    {
	      recEv->location = eventPos + nonAncestralLength - lastLength;
	      recEv->chrom = currChrom;
	      foundPosition=1;
	    }
	  lastPosition = tmp_anc->position;
	  tmp_anc = tmp_anc->next;
	}
      currChrom = currChrom->next;
    }
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
  assert(sumAnc != 0);
  chrtmp->next = newLeft;
  chrtmp = chrtmp->next;
  *noChrom = *noChrom + 1;
  
  currAnc = newRight->anc;
  sumAnc=0;
  while( currAnc != NULL)
    {
      sumAnc += currAnc->abits;
      currAnc = currAnc->next;
    }
  assert(sumAnc != 0);
  
  chrtmp->next = newRight;
  *noChrom = *noChrom + 1;
  *noChrom = *noChrom - 1;
}

chromosome* mergeChr(chromosome* ptrchr1, chromosome* ptrchr2)
{
  double epsilon = 1e-10;
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

void combineIdentAdjAncSegs(chromosome *ptrchr)
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

void coalescence(gsl_rng * r, unsigned int* noChrom, chrsample* chrom)
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

unsigned long long int ipow( unsigned long long int base, int exp)
{
  unsigned long long int result = 1;
  while( exp )
    {
      if ( exp & 1 )
        {
	  result *= (unsigned long long int)base;
        }
      exp >>= 1;
      base *= base;
    }
  return result;
}

int TestMRCAForAll(chrsample* chrom, unsigned int mrca)
{
  chromosome* currChrom = chrom->chrHead;
  ancestry* tmp_anc;
  while(currChrom != NULL)
    {
      tmp_anc = currChrom->anc;
      while(tmp_anc != NULL)
	{
	  if((tmp_anc->abits > 0)&&(tmp_anc->abits != mrca))
	    return(0);
	  tmp_anc = tmp_anc->next;
	}
      currChrom = currChrom->next;
    }
  return(1);
}

chrsample* create_sample(int noChrom)
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
