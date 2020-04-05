#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

gsl_rng * r;

/* options */
int opt_help = 0; /* print help message */
int opt_default = 0; /* print help message */

char version[] = "coalsim";

static void print_msg()
{
  printf("Usage: coalsim [OPTION]... \n Try 'coalsim -h' for more information.\n");
}

static void print_help()
{
  printf("Coalescent simulations: \n"
	 "-c <sample size>\n"
	 "-N <population size>\n"
	 "-r <recombination rate>\n"
	 "-s <seed for RNG>\n"
	 "-m <mutation rate>\n");
}

int main(int argc, char **argv)
{
  fillheader(version);
  show_header();
  double popSize = 10000;
  double recRate = 0.1;
  double mutRate = 3.0;
  unsigned int noSamples=5;
  unsigned int RGSeed=0;
  char* endPtr;
  opterr = 0;
  int c;
  const gsl_rng_type * T;

  while((c = getopt(argc, argv, "c:N:r:m:s:h")) != -1)
    switch(c)
      {
      case 'c':
	noSamples = strtoul(optarg,&endPtr,10);
	break;
      case 'N':
	popSize = atof(optarg);
	break;
      case 'r':
	recRate = atof(optarg);
	break;
      case 'm':
	mutRate = atof(optarg);
	break;
      case 's':
	RGSeed = strtoul(optarg,&endPtr,10);
	break;
      case 'h':
	opt_help = 1;
	break;
      case '?':
        if (isprint (optopt))
          fprintf(stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf(stderr,"Unknown option character `\\x%x'.\n",optopt);
        return 1;
      default:
	abort();
      }
  if(optind == 1)
    {
      print_msg();
      return 1;
    }
  if(opt_help)
    {
      print_help();
      return 1;
    }

 /* create a generator chosen by the
    environment variable GSL_RNG_TYPE */
  gsl_rng_env_setup();

  // gsl_rng_default_seed = 45567; 
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  if(RGSeed != 0)
    gsl_rng_set(r,RGSeed);
  
  unsigned int noChrom = noSamples;
  recombination_event recombEvent;
  chromosome* currentChrom = NULL;
  chrsample* chromSample = create_sample(noChrom);
  ancestry* tmp = NULL;
  int currChr=0;
  double ancLength=0;
  double eventLocation=0;
  double totalTime=0;
  double interArrivalTime=0;
  currentChrom = chromSample->chrHead;
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
      printf("interArrivaltime: %lf totalTime: %lf NoMut: %d NoRec: %d NoCoal: %d\n",interArrivalTime,totalTime,noMutations,noRec,noCoal);
    
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
  } 
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
  delete_sample(chromSample->chrHead);
  free(chromSample); 
}

