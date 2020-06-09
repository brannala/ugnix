#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

gsl_rng * r;

/* options */
int opt_help = 0; /* print help message */
int opt_default = 0; /* print help message */
int prn_chrom = 0; /* print chromosomes */
int prn_mrca = 0; /* print MRCAs */
int prn_mutations = 0; /* print mutations */
int prn_regions = 0; /* print all mrca regions */
int calc_mrca = 0; /* do mrca calculations */
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
	 "-u <specify scaling for bases: Mb Kb b>\n"
	 "-a <output mrca information: r=regions i=intervals s=summary>\n"
	 "-l <print detailed information about mutations>\n"
	 "-d <print chromosomes>\n"
	 "-m <mutation rate>\n");
}

int main(int argc, char **argv)
{
  fillheader(version);
  show_header();
  double popSize = 1000;
  double recRate = 0.05;
  double mutRate = 0.5;
  unsigned int noSamples=4;
  unsigned int RGSeed=0;
  char* endPtr;
  opterr = 0;
  int c;
  int seqUnits = 0;
  double cMtoMb = 1.0;
  const gsl_rng_type * T;

  while((c = getopt(argc, argv, "c:N:r:m:s:u:a:dlh")) != -1)
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
      case 'u':
	if(!strcmp("b",optarg))
	  seqUnits = 1;
	else
	  if((!strcmp("Kb",optarg))||(!strcmp("kb",optarg)))
	    seqUnits = 2;
	  else
	    if((!strcmp("Mb",optarg))||(!strcmp("mb",optarg)))
	      seqUnits = 3;
	    else
	      {
		fprintf(stderr,"Unknown specifier '%s' for -u.\n",optarg);
		return 1;
	      }
	break;
      case 'a':
	if(!strcmp("r",optarg))
	  {
	    prn_regions = 1;
	    calc_mrca = 1;
	  }
	else
	  if(!strcmp("i",optarg))
	    {
	    prn_mrca = 1;
	    calc_mrca = 1;
	    }
	  else
	    if(!strcmp("s",optarg))
	      {
		calc_mrca = 1;
	      }
	    else
	      {
		fprintf(stderr,"Unknown specifier '%s' for -a.\n",optarg);
		return 1;
	      }
	break;
      case 'd':
	prn_chrom = 1;
	break;
      case 'l':
	prn_mutations = 1;
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
  mutation* mutation_list = NULL;
  double ancLength=0;
  double eventLocation=0;
  double totalTime=0;
  double interArrivalTime=0;
  currentChrom = chromSample->chrHead;
  double totRate=0;
  double coalProb = 0;
  double recProb = 0;
  double smalldiff = 1e-8;
  int noMutations=0;
  int noRec=0;
  int noCoal=0;
  struct mrca_list* head;
  struct mrca_summary* mrca_head = NULL;
  unsigned int mrca=0;
  for(unsigned int i=0; i <noSamples; i++)
    {
      mrca += ipow(2,i);
    }
  long chromTotBases = recRate*100*cMtoMb*1e6;
  char* baseUnit = malloc(sizeof(char)*5);
  if(seqUnits == 1)
    strcpy(baseUnit,"bps");
  else
    if(seqUnits == 2)
      strcpy(baseUnit,"kb");
    else
      if(seqUnits == 3)
	strcpy(baseUnit,"Mb");
      else
	strcpy(baseUnit,"");

  /* main simulation loop */
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
    totalTime += interArrivalTime;
    prob = gsl_rng_uniform_pos(r);
    if(prob <= coalProb)
      /* coalescence event */
      {

	coalescent_pair pair;
	getCoalPair(r,noChrom,&pair);
	coalescence(pair,&noChrom, chromSample);
      	if(calc_mrca)
	  getMRCAs(&head,currentChrom,chromSample,totalTime,mrca);
	noCoal++;
      }
    else
      if(prob <= (coalProb + recProb))
	/* recombination event */
	{
	  eventLocation = ancLength*gsl_rng_uniform_pos(r);
	  assert(eventLocation <= ancLength);
	  getRecEvent(chromSample, eventLocation, &recombEvent);
	  recombination(&noChrom,recombEvent,chromSample);
	  noRec++;
	}
      else
	/* mutation event */
	{
	  mutation* tmpMut = malloc(sizeof(mutation));
	  mutation* mcurr;
	  tmpMut->next = NULL;
	  eventLocation = ancLength*gsl_rng_uniform_pos(r);
	  assert(eventLocation <= ancLength);
	  getMutEvent(chromSample, eventLocation, tmpMut, totalTime);
	  if(mutation_list == NULL)
	    mutation_list = tmpMut;
	  else
	    {
	      mcurr = mutation_list;
	      while(mcurr->next != NULL)
		mcurr = mcurr->next;
	      mcurr->next = tmpMut;
	    }
	  noMutations++;
	}
  }

  /* summarize run input and output */
  printf("N:%.0f n:%d r:%.2f ",popSize,noSamples,recRate);
  printf("Mutation_Rate: %.2f/Chr",mutRate);
  if(seqUnits)
    printf(" %.3e/base",mutRate/chromTotBases);
  printf("\n");
  if(seqUnits)
    printf("Chr_Length: %ld%s (Assumes %.2fcM/Mb)\n",
	   convertToBases(chromTotBases,seqUnits,1),baseUnit,cMtoMb);
  printf("No_Recombinations: %d ",noRec);
  printf("No_Mutations: %d ",noMutations);
  printf("No_Ancestral_Chromosomes: %d\n",noChrom);
  printf("Oldest_TMRCA: %.2lf ",totalTime);
  
  if(calc_mrca)
    MRCAStats(head,mrca_head,smalldiff,chromTotBases,seqUnits,baseUnit,prn_mrca,prn_regions);
  else
    printf("\n");

  if(prn_mutations)
    {
      printMutations(mutation_list,chromTotBases,seqUnits,baseUnit,noSamples,mrca);
    }
  
  if(prn_chrom)
    {
      printChromosomes(chromSample,noSamples);
    }
  
  delete_sample(chromSample->chrHead);
  free(chromSample); 
}

