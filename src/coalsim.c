#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

FILE *out_file;

gsl_rng * r;

/* options */
int opt_help = 0; /* print help message */
int opt_default = 0; /* print help message */
int prn_chrom = 0; /* print chromosomes */
int prn_mrca = 0; /* print MRCAs */
int prn_mutations = 0; /* print mutations */
int prn_regions = 0; /* print all mrca regions */
int calc_mrca = 0; /* do mrca calculations */
int prn_sequences = 0; /* print sequence data to output file */
int prn_genetrees = 0; /* print gene trees for mrca regions */
int gtrees_to_stdout = 0; /*print gene trees to stdout */
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
	 "-g <print gene trees for mrca regions s=screen f=file>\n"
	 "-o <sequence output file name>\n"
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
  char* outfile = malloc(sizeof(char)*30);
  unsigned int noSamples=4;
  unsigned int RGSeed=0;
  char* endPtr;
  opterr = 0;
  int c;
  int seqUnits = 0;
  double cMtoMb = 1.0;
  const gsl_rng_type * T;
  

  while((c = getopt(argc, argv, "c:N:r:m:s:u:a:o:g:dlh")) != -1)
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
      case 'o':
	prn_sequences = 1;
	strncpy(outfile,optarg,30);
	out_file = fopen(outfile,"w");
	if(out_file == NULL)
	  {
	    fprintf(stderr, "Error: failed to open output file %s\n",outfile);
	    exit(1);
	  }
	break;
      case 'g':
	prn_genetrees = 1;
	if(!strcmp("s",optarg))
	   gtrees_to_stdout = 1;
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
  chrsample* chromSample = create_sample(noChrom);
  struct coalescent_events* coalescent_list = NULL;
  mutation* mutation_list = NULL;
  char** sequences;
  double ancLength=0;
  double eventLocation=0;
  double totalTime=0;
  double interArrivalTime=0;
  double totRate=0;
  double coalProb = 0;
  double recProb = 0;
  double smalldiff = 1e-8;
  int noMutations=0;
  int noRec=0;
  int noCoal=0;
  struct mrca_list* head = NULL;
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
	if(prn_genetrees)
	    updateCoalescentEvents(&coalescent_list,chromSample,totalTime);
      	if(calc_mrca || prn_genetrees) 
	  getMRCAs(&head,chromSample,totalTime,mrca);
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
  printf("Mutation_Rate: %.3e/Chr",mutRate);
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

  if(prn_sequences)
    {
      sequences = simulateSequences(mutation_list,chromTotBases,noSamples,r);
      for(int i=0; i<noSamples; i++)
	{
	  fprintf(out_file,">sample%d\n",i);
	  for(int j=0; j<chromTotBases; j++)
	    fprintf(out_file,"%c",sequences[i][j]);
	  fprintf(out_file,"\n");
	}
    }

  if(prn_genetrees)
    {
      struct mrca_list* tmp_mrca_list = head;
      struct tree* t1;
      struct geneTree* g2;
      printf("\nGene trees for chromosome regions\n");
      printf("-----------------------------------\n");	
      printf("MRCA Interval:       Gene Tree:\n\n");
      while(tmp_mrca_list != NULL)
	{
	  struct geneTree* g1 = getGeneTree(tmp_mrca_list->lower_end, tmp_mrca_list->upper_end,coalescent_list);
	  g2=g1;
	  t1 = malloc(sizeof(struct tree));
	  t1->abits = g2->abits;
	  t1->time = g2->time;
	  t1->left = NULL;
	  t1->right = NULL;
	  g2 = g2->next;
	  while(g2 != NULL)
	    {
	      struct tree* tmp = t1;
	      addNode(g2->abits,g2->time,tmp);
	      g2 = g2->next;
	    }
	  fillTips(t1);
	  printf("(%f, ", tmp_mrca_list->lower_end);
	  printf("%f) ", tmp_mrca_list->upper_end);
	  printTree(t1,noSamples);
	  printf("\n");
	  g2=g1;
	  g1 = g1->next;
	  while(g1 != NULL)
	    {
	      free(g2);
	      g2=g1;
	      g1 = g1->next;
	    }
	  free(g2);
	  tmp_mrca_list = tmp_mrca_list->next;
	}
    }
  
  if(prn_chrom)
    {
      printChromosomes(chromSample,noSamples);
    }

  /* clean up memory */
  delete_sample(chromSample->chrHead);
  free(chromSample); 
  free(baseUnit);
  free(outfile);
  mutation* mcurr = mutation_list;
  while(mutation_list != NULL)
    {
      mcurr = mutation_list;
      mutation_list = mutation_list->next;
      free(mcurr);
    }
  struct mrca_list* hcurr = head;
  if(calc_mrca)
    {
      while(head != NULL)
	{
	  hcurr = head;
	  head = head->next;
	  free(hcurr);
	}
      struct mrca_summary* mhcurr = mrca_head;
      while(mrca_head != NULL)
	{
	  mhcurr = mrca_head;
	  mrca_head = mrca_head->next;
	  free(mhcurr);
	}
    }

  if(out_file != NULL)
    fclose(out_file);
}

