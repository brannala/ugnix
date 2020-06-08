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
	 "-d print chromosomes\n"
	 "-a print MRCAs\n"
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
  ancestry* tmp = NULL;
  mutation* mutation_list = NULL;
  int currChr=0;
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
  struct mrca_list* curr;
  double newlower=0;
  double newupper=0;
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
    /*    if(!(eventNo % 5000))
	  printf("interArrivaltime: %lf totalTime: %lf NoMut: %d NoRec: %d NoCoal: %d\n",interArrivalTime,totalTime,noMutations,noRec,noCoal); */
    
    totalTime += interArrivalTime;
    prob = gsl_rng_uniform_pos(r);
    if(prob <= coalProb)
      /* coalescence event */
      {

	coalescent_pair pair;
	getCoalPair(r,noChrom,&pair);
	coalescence(pair,&noChrom, chromSample);

      	if(calc_mrca)
	  {
	    // collect information on intervals and ages of unique mrca's
	    int firstInt=1;
	    currentChrom = chromSample->chrHead; 
	    while(currentChrom->next != NULL)
	      currentChrom = currentChrom->next;
	    tmp = currentChrom->anc;
	    while((tmp->next != NULL)||firstInt)
	      {
		if(firstInt)
		  {
		    if(tmp->abits == mrca)
		      {
			newlower = 0.0;
			newupper = tmp->position;
			addMRCAInterval(&head,newlower,newupper,totalTime);
		      }
		    firstInt=0;
		    if(tmp->next != NULL)
		      tmp = tmp->next;
		  }
		else
		  {
		    if(tmp->next->abits == mrca)
		      {
			newlower = tmp->position;
			newupper = tmp->next->position;
			addMRCAInterval(&head,newlower,newupper,totalTime);
		      }
		    tmp = tmp->next;
		  }
	      } 
	  } 
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
  printf("N:%.0f n:%d r:%.2f ",popSize,noSamples,recRate);
  printf("Mutation_Rate: %.2f/Chr",mutRate);
  if(seqUnits)
    printf(" %.3e/base",mutRate/chromTotBases);
  printf("\n");
  printf("Chr_Length: %ld%s (Assumes %.2fcM/Mb)\n",
  	 convertToBases(chromTotBases,seqUnits,1),baseUnit,cMtoMb);
  printf("No_Recombinations: %d ",noRec);
  printf("No_Mutations: %d ",noMutations);
  printf("No_Ancestral_Chromosomes: %d\n",noChrom);
  printf("Oldest_TMRCA: %.2lf ",totalTime);
  
  if(prn_mutations)
    {
      mutation* tmp1 = mutation_list;
      while(tmp1 != NULL)
	{
	  if(tmp1->abits != mrca)
	    {
	      if(seqUnits)
		printf("Location: %ld%s Age: %f Chrom: ",
		       convertToBases(chromTotBases,seqUnits,tmp1->location),baseUnit,tmp1->age);
	      else
		printf("Location: %f Age: %f Chrom: ",
		       tmp1->location,tmp1->age);
	      displayBits(tmp1->abits,noSamples);
	      printf("\n");
	    }
	  tmp1 = tmp1->next;	
	}
    }
  
  if(calc_mrca)
    {

      curr=head;
      int firstEntry=1;
      while(curr != NULL)
	{
	  struct mrca_summary* curr_sum;
	  struct mrca_summary* last_sum;
	  struct mrca_summary* tmp;
	  int isHead=1;
	  curr_sum = mrca_head;
	  while((curr_sum != NULL)&&(curr_sum->age < curr->age))
	    {
	      last_sum = curr_sum;
	      curr_sum = curr_sum->next;
	      isHead = 0;
	    }
	  if(curr_sum == NULL)
	    {
	      if(firstEntry)
		{
		  curr_sum = malloc(sizeof(struct mrca_summary));
		  curr_sum->age = curr->age;
		  curr_sum->length = curr->upper_end - curr->lower_end;
		  curr_sum->numInts = 1;
		  curr_sum->next = NULL;
		  mrca_head = curr_sum;
		  firstEntry = 0;
		  curr = curr->next;
		  continue;
		}
	      else
		{
		  last_sum->next = malloc(sizeof(struct mrca_summary));
		  last_sum->next->age = curr->age;
		  last_sum->next->length = curr->upper_end - curr->lower_end;
		  last_sum->next->numInts = 1;
		  last_sum->next->next = NULL;
		  curr = curr->next;
		  continue;
		}
	    }
	  if(fabs(curr_sum->age - curr->age) < smalldiff)
	    {
	      curr_sum->numInts++;
	      curr_sum->length += curr->upper_end - curr->lower_end;
	      curr = curr->next;
	      continue;
	    }
	  else
	    {
	      tmp = malloc(sizeof(struct mrca_summary));
	      tmp->age = curr->age;
	      tmp->length = curr->upper_end - curr->lower_end;
	      tmp->numInts = 1;
	      tmp->next = curr_sum;
	      if(isHead)
		mrca_head = tmp;
	      else
		last_sum->next = tmp;
	    }
	   curr = curr->next;
	}

      struct mrca_summary* curr_sum;
      double mean_tmrca=0;
      double sumsqr=0;
      double v2=0;
      int no_segs=0;
      double largest_mrca=0;
      double smallest_mrca=1e9;
      double youngest_mrca=1e9;
      curr_sum = mrca_head;
      /* calculate mean and var of tmrca across segments */
      while(curr_sum != NULL)
	{
	  no_segs++;
	  mean_tmrca += curr_sum->age;
	  if(largest_mrca < curr_sum->length)
	    largest_mrca = curr_sum->length;
	  if(smallest_mrca > curr_sum->length)
	    smallest_mrca = curr_sum->length;
	  if(youngest_mrca > curr_sum->age)
	     youngest_mrca = curr_sum->age;
	  curr_sum = curr_sum->next;
	}
      curr_sum = mrca_head;
      no_segs=0;
      while(curr_sum != NULL)
	{
	  no_segs++;
	  sumsqr += (curr_sum->age - mean_tmrca)*(curr_sum->age - mean_tmrca);
	  v2 = (1.0/no_segs)*(curr_sum->age - mean_tmrca)*(curr_sum->age - mean_tmrca) + ((no_segs-1.0)/no_segs)*v2;
	  curr_sum = curr_sum->next;
	}
      printf("Youngest_TMRCA: %.2f ",youngest_mrca);
      printf("Avg_TMRCA: %.2f\n",mean_tmrca/no_segs);
      printf("No_MRCAs: %d ",no_segs);

      if(seqUnits)
	printf("Largest_MRCA: %ld%s ",
	       convertToBases(chromTotBases,seqUnits,largest_mrca),baseUnit);
      else
      printf("Largest_MRCA: %.2f ",largest_mrca);
      if(seqUnits)
	printf("Smallest_MRCA: %ld%s.\n",
	       convertToBases(chromTotBases,seqUnits,smallest_mrca),baseUnit);
      else
      printf("Smallest_MRCA: %.2f.\n",smallest_mrca);



      // printf("SE[Mean(TMRCA)]: %f",sqrt(sumsqr/no_segs)/sqrt(no_segs*1.0));
      // printf("SE[V2]: %f\n\n",sqrt(v2)/sqrt(no_segs*1.0));
      if(prn_mrca)
	{
	  printf("TMRCA Summary\n");
	  printf("-----------------------------\n\n");	
	  curr_sum = mrca_head;
	  while(curr_sum != NULL)
	    {
	      if(seqUnits)
		printf("Length: %ld%s tmrca: %f\n",
		       convertToBases(chromTotBases,seqUnits,curr_sum->length),baseUnit,curr_sum->age);
	      else
		printf("Length: %f tmrca: %f\n, ", curr_sum->length, curr_sum->age);
	      curr_sum = curr_sum->next;
	    }
	}
      if(prn_regions)
	{
	  printf("\nTMRCAs for chromosome regions\n");
	  printf("-----------------------------\n\n");	
	  curr=head;
	  while(curr != NULL)
	    {
	      printf("(%f, ", curr->lower_end);
	      printf("%f) ", curr->upper_end);
	      printf(" tmrca: %f\n",curr->age);
	      curr = curr->next;
	    }
	}
    }
    
  if(prn_chrom)
    {
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
    }
  delete_sample(chromSample->chrHead);
  free(chromSample); 
}

