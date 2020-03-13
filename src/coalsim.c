#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

gsl_rng * r;


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
  chromosome* currentChrom = NULL;
  chrsample* chromSample = create_sample(noChrom);
  ancestry* tmp = NULL;
  int currChr=0;
  double ancLength=0;
  double eventLocation=0;
  double totalTime=0;
  double interArrivalTime=0;
  currentChrom = chromSample->chrHead;
  double popSize = 10000;
  double recRate = 0.1;
  double mutRate = 3.0;
  double totRate=0;
  double coalProb = 0;
  double recProb = 0;
  
  int noMutations=0;
  int noRec=0;
  int noCoal=0;
  unsigned int mrca=0;
  for(unsigned int i=0; i <noSamples; i++)
    {
      mrca += ipow(2,i);
    }


  int eventNo = 0;
  while((noChrom > 1) && (!TestMRCAForAll(chromSample, mrca)) )
  {
    double prob = 0;
    eventNo++;
    ancLength = totalAncLength(chromSample);
    assert(ancLength <= noSamples);
    totRate = (noChrom*(noChrom-1)/2.0)*(1.0/(2.0*popSize))+(recRate +mutRate)*ancLength;
    coalProb = ((noChrom*(noChrom-1)/2.0)*(1.0/(2.0*popSize)))/totRate;
    recProb = recRate*ancLength/totRate;
    assert(coalProb + recProb < 1.0);
    interArrivalTime = gsl_ran_exponential(r, 1.0/totRate);
    if(!(eventNo % 5000))
      printf("interArrivaltime: %lf NoRec: %d NoCoal: %d\n",interArrivalTime,noRec,noCoal);
    
    totalTime += interArrivalTime;
    prob = gsl_rng_uniform_pos(r);
    if(prob <= coalProb)
      {
	coalescent_pair pair;
	getCoalPair(r,noChrom,&pair);
	coalescence(pair,&noChrom, chromSample);
	noCoal++;
      }
    else
      if(prob <= (coalProb + recProb))
	{
	  eventLocation = ancLength*gsl_rng_uniform_pos(r);
	  assert(eventLocation <= ancLength);
	  getRecEvent(chromSample, eventLocation, &recombEvent);
	  recombination(&noChrom,recombEvent,chromSample);
	  noRec++;
	}
      else
	noMutations++;

    /*    printf("recNo: %d\n",noRec);
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
	} */
    
  } 
  /*  ancLength = totalAncLength(chromSample);
  eventLocation = ancLength*gsl_rng_uniform_pos(r);
  getRecEvent(chromSample, eventLocation, &recombEvent);
  recombination(&noChrom,recombEvent,chromSample); */
  printf("recNo: %d\n",noRec);
  printf("mutNo: %d\n",noMutations);
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



  
  // coalescence(&noChrom, chromSample);
  /* currChr=0;
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
  printf("\nNoChrom: %d\n",noChrom); */
  //  printf("U(0,1): %lf\n",gsl_rng_uniform_pos(r));
  // printf("U(0,1): %lf",gsl_ran_flat(r,0,1.0));
  delete_sample(chromSample->chrHead);
  free(chromSample); 
}

