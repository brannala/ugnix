#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

struct geneTree* SortedMerge(struct geneTree* a, struct geneTree* b) 
{ 
	struct geneTree* result = NULL; 

	if (a == NULL) 
		return (b); 
	else if (b == NULL) 
		return (a); 
	if (a->time >= b->time) { 
		result = a; 
		result->next = SortedMerge(a->next, b); 
	} 
	else { 
		result = b; 
		result->next = SortedMerge(a, b->next); 
	} 
	return (result); 
} 

void FrontBackSplit(struct geneTree* source, struct geneTree** frontRef, struct geneTree** backRef) 
{ 
	struct geneTree* fast; 
	struct geneTree* slow; 
	slow = source; 
	fast = source->next; 
	while (fast != NULL) { 
		fast = fast->next; 
		if (fast != NULL) { 
			slow = slow->next; 
			fast = fast->next; 
		} 
	} 
	*frontRef = source; 
	*backRef = slow->next; 
	slow->next = NULL; 
} 

/* sorts the linked list by changing next pointers (not data) */
void MergeSort(struct geneTree** headRef) 
{ 
	struct geneTree* head = *headRef; 
	struct geneTree* a; 
	struct geneTree* b; 


	if ((head == NULL) || (head->next == NULL)) { 
		return; 
	} 

	FrontBackSplit(head, &a, &b); 
	MergeSort(&a); 
	MergeSort(&b); 
	*headRef = SortedMerge(a, b); 
} 

int isSingleton(unsigned int x)
{
  if(x == 1)
    return 1;
  else
    {
      unsigned int z=1;
      while(z < 1073741825)
	if(x == (z*=2))
	  return 1;
    }
  return 0;
}

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

chromosome* copy_chrom(chromosome* sourceChr)
{
  chromosome* newChr;
  ancestry* currNew;
  ancestry* currOld;
  newChr = malloc(sizeof(chromosome));
  newChr->next = NULL;
  currOld = sourceChr->anc;
  newChr->anc = malloc(sizeof(ancestry));
  int firstAnc=1;
  currNew = newChr->anc;
  while(currOld != NULL)
    {
      if(!firstAnc)
	{
	  currNew->next = malloc(sizeof(ancestry));
	  currNew = currNew->next;
	}
      firstAnc = 0;
      currNew->position = currOld->position;
      currNew->abits = currOld->abits;
      currOld = currOld->next;
      
    }
  currNew->next = NULL;
  return newChr;
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
  /* Net chromosome count change: -1 (deleted original) + 2 (newLeft + newRight) = +1
     The increment at line 305 after adding newLeft accounts for this */
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
      else
	if(tmp->next != NULL)
	  tmp = tmp->next;
    }
}

void getCoalPair(gsl_rng * r, unsigned int noChrom, coalescent_pair* pair)
{
  if(noChrom > 2)
    {
      /* Pick first chromosome uniformly from [0, noChrom-1] */
      pair->chr1 = gsl_rng_uniform_int(r, noChrom);
      /* Pick second chromosome uniformly from remaining [0, noChrom-2] */
      pair->chr2 = gsl_rng_uniform_int(r, noChrom - 1);
      /* Adjust to skip chr1 */
      if(pair->chr2 >= pair->chr1)
	pair->chr2++;
    }
  else
    {
      pair->chr1 = 0;
      pair->chr2 = 1;
    }
}

void coalescence(coalescent_pair pair, unsigned int* noChrom, chrsample* chrom)
{

  *noChrom = *noChrom - 1;
  chromosome* tmp;
  chromosome* ptrchr1 = NULL;
  chromosome* ptrchr2 = NULL;
  chromosome* commonAnc = NULL;
  ptrchr1 = getChrPtr(pair.chr1, chrom);
  ptrchr2 = getChrPtr(pair.chr2, chrom);
  commonAnc = mergeChr(ptrchr1, ptrchr2);
  // check?
  // combineIdentAdjAncSegs(commonAnc);
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

void updateCoalescentEvents(struct coalescent_events** coalescent_list, chrsample* chromSample, double totalTime)
{
	    struct coalescent_events* tmpCList;
	    chromosome* tmpc;
	    if(*coalescent_list == NULL)
	      {
		*coalescent_list = malloc(sizeof(struct coalescent_events));
		tmpCList = *coalescent_list;
	      }
	    else
	      {
		tmpCList = *coalescent_list;
		while(tmpCList->next != NULL)
		  tmpCList = tmpCList->next;
		tmpCList->next = malloc(sizeof(struct coalescent_events));
		tmpCList = tmpCList->next;
	      }

	    tmpc = chromSample->chrHead;
	    while(tmpc->next != NULL)
	      tmpc = tmpc->next;
	    tmpCList->chr = copy_chrom(tmpc);
	    combineIdentAdjAncSegs(tmpCList->chr);
	    tmpCList->time = totalTime;
	    tmpCList->next = NULL;
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

struct geneTree* getGeneTree(double lower, double upper, struct coalescent_events* coalescent_list, unsigned int mrca)
{
  struct geneTree* geneT = NULL;
  struct geneTree* currGT = NULL;
  struct coalescent_events* localCL=coalescent_list;
  ancestry* localAnc;

  #ifdef DEBUG_GENETREE
  fprintf(stderr, "DEBUG getGeneTree: lower=%.4f upper=%.4f mrca=%u\n", lower, upper, mrca);
  struct coalescent_events* debugCL = coalescent_list;
  int eventNum = 0;
  while(debugCL != NULL) {
    fprintf(stderr, "  Event %d: time=%.2f, chr->anc->abits=%u\n",
            eventNum++, debugCL->time, debugCL->chr->anc->abits);
    debugCL = debugCL->next;
  }
  #endif

  while((localCL != NULL)&&((currGT == NULL)||(currGT->abits != mrca)))
    {
      localAnc = localCL->chr->anc;
      double lastPos=0.0;
      int foundAnc=0;
      while((localAnc != NULL) && !foundAnc)
	{
	  if(!(isSingleton(localAnc->abits) || localAnc->abits == 0))
	    {
	      if((upper <= localAnc->position)&&(lower >= lastPos))
		{
		  struct geneTree* tmpGT = geneT;
		  while((tmpGT != NULL) && !foundAnc)
		    {
		      if(tmpGT->abits == localAnc->abits)
			foundAnc=1;
		      tmpGT = tmpGT->next;
		    }
		  if(!foundAnc)
		    {
		      if(geneT == NULL)
			{
			  geneT = malloc(sizeof(struct geneTree));
			  geneT->next = NULL;
			  currGT = geneT;
			}
		      else
			{
			  currGT->next = malloc(sizeof(struct geneTree));
			  currGT = currGT->next;
			  currGT->next = NULL;
			}
		      currGT->abits = localAnc->abits;
		      currGT->time = localCL->time;
		      foundAnc = 1;
		    }
		}
	    }
	  lastPos = localAnc->position;
	  localAnc = localAnc->next;
	}
      localCL = localCL->next;
    }
  MergeSort(&geneT);
  return geneT;
}

void addNode(unsigned int val, double time, struct tree* lroot)
{
  while(lroot->right != NULL)
    {
      if(lroot->right->abits == val)
	{
	  lroot->right->time = time;
	  return;
	}
      else
	if(lroot->left->abits == val)
	{
	  lroot->left->time = time;
	  return;
	}  
      if((lroot->right->abits & val)!=0)
	lroot = lroot->right;
      else
	lroot = lroot->left;
    }
  lroot->right = malloc(sizeof(struct tree));
  lroot->right->right = NULL;
  lroot->right->left = NULL;
  lroot->right->abits = val;
  lroot->right->time = time;
  lroot->left = malloc(sizeof(struct tree));
  lroot->left->abits = ~val & lroot->abits;
  lroot->left->left = NULL;
  lroot->left->right = NULL;
}

void splitNode(struct tree* lroot)
{
  unsigned int bitmask = 1;
  while(!(bitmask & lroot->abits))
    bitmask = bitmask << 1;
  if((~bitmask & lroot->abits)!=0)
    {
      lroot->left = malloc(sizeof(struct tree));
      lroot->right = malloc(sizeof(struct tree));
      lroot->left->left = NULL;
      lroot->left->right = NULL;
      lroot->left->abits = bitmask;
      lroot->right->right = NULL;
      lroot->right->left = NULL;
      lroot->right->abits = (~bitmask & lroot->abits);
    }
  return;
}

void fillTips(struct tree* lroot)
{
  if(lroot->left != NULL)
    {
      fillTips(lroot->left);
      fillTips(lroot->right);
    }
  else
    splitNode(lroot);
  return;
}

unsigned int binaryToChrLabel(unsigned int x, int noSamples)
{
  unsigned int bitmask = 1;
  unsigned int pos = 1;
  while(!(bitmask & x) && (pos <= noSamples))
    {
      bitmask = bitmask << 1;
      pos++;
    }
  return pos;

}

/* Helper function to print tree recursively with proper branch lengths */
static void printTreeNewickHelper(struct tree* node, int noSamples, int toScreen,
                                   FILE* tree_file, double parentTime)
{
  FILE* out = toScreen ? stderr : tree_file;

  if(node->left != NULL && node->right != NULL)
    {
      /* Internal node with two children */
      fprintf(out, "(");
      printTreeNewickHelper(node->left, noSamples, toScreen, tree_file, node->time);
      fprintf(out, ",");
      printTreeNewickHelper(node->right, noSamples, toScreen, tree_file, node->time);
      fprintf(out, ")");
      /* Branch length from this node to parent */
      double branchLen = parentTime - node->time;
      if(branchLen > 0)
        fprintf(out, ":%.4f", branchLen);
    }
  else
    {
      /* Tip node - branch length is from time 0 to parent */
      fprintf(out, "%d:%.4f", binaryToChrLabel(node->abits, noSamples), parentTime);
    }
}

/* Print gene tree in Newick format with branch lengths in generations */
void printTreeNewick(struct tree* root, int noSamples, int toScreen, FILE* tree_file)
{
  FILE* out = toScreen ? stderr : tree_file;

  if(root == NULL) return;

  if(root->left != NULL && root->right != NULL)
    {
      fprintf(out, "(");
      printTreeNewickHelper(root->left, noSamples, toScreen, tree_file, root->time);
      fprintf(out, ",");
      printTreeNewickHelper(root->right, noSamples, toScreen, tree_file, root->time);
      fprintf(out, ");");  /* Root has no branch length, ends with semicolon */
    }
  else
    {
      /* Degenerate case: single tip */
      fprintf(out, "%d;", binaryToChrLabel(root->abits, noSamples));
    }
}

/* Free tree memory recursively */
void freeTree(struct tree* node)
{
  if(node == NULL) return;
  freeTree(node->left);
  freeTree(node->right);
  free(node);
}

/* Legacy function - kept for compatibility but uses absolute times (deprecated) */
void printTree(struct tree* lroot, int noSamples, int toScreen, FILE* tree_file)
{
  if(lroot->left != NULL)
    {
      toScreen? fprintf(stderr,"(") : fprintf(tree_file,"(");
      printTree(lroot->left,noSamples,toScreen,tree_file);
      if(lroot->left->left == NULL)
	toScreen? fprintf(stderr,",") : fprintf(tree_file,",");
    }
  if(lroot->right != NULL)
    {
      printTree(lroot->right,noSamples,toScreen,tree_file);
      toScreen? fprintf(stderr,")") : fprintf(tree_file,")");
      toScreen? fprintf(stderr,":%.2f",lroot->time) : fprintf(tree_file,":%.2f",lroot->time);
    }
  if(lroot->left == NULL)
    toScreen? fprintf(stderr,"%d",binaryToChrLabel(lroot->abits,noSamples)) : fprintf(tree_file,"%d",binaryToChrLabel(lroot->abits,noSamples));
  return;
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
	    return(0); // not all zeros or all ones therefore not mrca of sample
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

void addMRCAInterval(struct mrca_list** head, double newlower,
				  double newupper, double newage)
{
  struct mrca_list* current_intv;
  struct mrca_list* last_intv=NULL;
  struct mrca_list* new_intv;
  int isHead = 1;
  // double smalldiff = 1e-8;
  double currNewlower = newlower;
  assert(newupper > newlower);
  current_intv = *head;
  while((current_intv != NULL)&&(current_intv->upper_end <= newlower) /* &&(current_intv->next != NULL) */)
    /* ! new           L-----U
       ! curr L-----U                   */
    /* found first overlapping interval */
     {
      last_intv = current_intv;
      current_intv = current_intv->next;
      isHead = 0;
    }
  
  while(1)
    {
      if((current_intv == NULL)||(newupper <= current_intv->lower_end))
	/* new L------U           OR        new   L-------U
	   curr         L------U      last L------U          */    
	{
	  new_intv = malloc(sizeof(struct mrca_list));
	  new_intv->age = newage;
	  new_intv->lower_end = currNewlower;
	  new_intv->upper_end = newupper;
	  new_intv->next = current_intv;
	  if(isHead)
	    *head = new_intv;
	  else
	    last_intv->next = new_intv;
	  return;
	}
      else
	if(currNewlower < current_intv->lower_end)
	  /* new L----->
             curr   L-----U */
	  {
	    new_intv = malloc(sizeof(struct mrca_list));
	    new_intv->age = newage;
	    new_intv->lower_end = currNewlower;
	    new_intv->upper_end = current_intv->lower_end;
	    new_intv->next = current_intv;
	    if(isHead)
	      *head = new_intv;
	    else
	      last_intv->next = new_intv;
	  }
      if(newupper <= current_intv->upper_end)
	return;
      else
	if(currNewlower <= current_intv->upper_end)
	  currNewlower = current_intv->upper_end;
      last_intv = current_intv;
      current_intv = current_intv->next;
    }
} 
 

void getMRCAs(struct mrca_list** head, chrsample* chromSample, double totalTime, unsigned int mrca)
{
  double newlower=0;
  double newupper=0;
  ancestry* tmp = NULL;
  chromosome* currentChrom;
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
	      addMRCAInterval(head,newlower,newupper,totalTime);
	    }
	  firstInt=0;
	  //  if(tmp->next != NULL)
	  //  tmp = tmp->next;
	}
      else
	{
	  if(tmp->next->abits == mrca)
	    {
	      newlower = tmp->position;
	      newupper = tmp->next->position;
	      addMRCAInterval(head,newlower,newupper,totalTime);
	    }
	  tmp = tmp->next;
	}
    } 
}

void MRCAStats(struct mrca_list* head, struct mrca_summary* mrca_head, double smalldiff, long chromTotBases,
		  int seqUnits, char* baseUnit, int prn_mrca, int prn_regions)
{
  struct mrca_list* curr;
  struct mrca_summary* curr_sum;
  struct mrca_summary* last_sum;
  struct mrca_summary* tmp;
  double mean_tmrca=0;
  int no_segs=0;
  double largest_mrca=0;
  double smallest_mrca=1e9;
  double youngest_mrca=1e9;
  
  curr=head;
  int firstEntry=1;
  while(curr != NULL)
    {
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
        
  /* calculate summary statistics of tmrca/mrca across segments */
  curr_sum = mrca_head;
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
  printf("Youngest_TMRCA: %.2f ",youngest_mrca);
  printf("Avg_TMRCA: %.2f\n",mean_tmrca/no_segs);
  printf("No_MRCAs: %d ",no_segs);
  if(seqUnits)
    printf("Largest_MRCA: %ld%s ",
	   convertToBases(chromTotBases,seqUnits,largest_mrca),baseUnit);
  else
    printf("Largest_MRCA: %.6f ",largest_mrca);
  if(seqUnits)
    printf("Smallest_MRCA: %ld%s.\n",
	   convertToBases(chromTotBases,seqUnits,smallest_mrca),baseUnit);
  else
    printf("Smallest_MRCA: %.6f.\n",smallest_mrca);

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
	    printf("Length: %f tmrca: %f\n", curr_sum->length, curr_sum->age);
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

void getMutEvent(chrsample* chrom, double eventPos, mutation* mutEv, double time)
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
	      mutEv->location = eventPos + nonAncestralLength - lastLength;
	      mutEv->age = time;
	      mutEv->abits = tmp_anc->abits;
	      foundPosition=1;
	    }
	  lastPosition = tmp_anc->position;
	  tmp_anc = tmp_anc->next;
	}
      currChrom = currChrom->next;
    }
}

void printMutations(mutation* mutation_list, long chromTotBases, int seqUnits,
		    char* baseUnit, unsigned int noSamples, unsigned int mrca)
{
  printf("\nMutation List\n------------------------------------\n");
  mutation* tmp1 = mutation_list;
  while(tmp1 != NULL)
    {
      if(tmp1->abits != mrca)
	{
	  if(seqUnits)
	    printf("Pos: %ld%s Age: %.2f Chrom: ",
		   convertToBases(chromTotBases,seqUnits,tmp1->location),baseUnit,tmp1->age);
	  else
	    printf("Pos: %f Age: %.2f Chrom: ",
		   tmp1->location,tmp1->age);
	  displayBits(tmp1->abits,noSamples);
	  printf("\n");
	}
      tmp1 = tmp1->next;	
    }
}

void printChromosomes(chrsample* chromSample, unsigned int noSamples)
{
  int currChr=0;
  chromosome* currentChrom=NULL;
  ancestry* tmp=NULL;
  currentChrom = chromSample->chrHead;
  while(currentChrom != NULL)
    {
      printf("\nChr: %d: \n",currChr);
      tmp = currentChrom->anc;
      while(tmp != NULL)
	{
	  if((tmp->next == NULL) || (tmp->next->abits != tmp->abits))
	    { 
	      displayBits(tmp->abits,noSamples);
	      printf(" %lf \n",tmp->position);  
	      tmp = tmp->next;
	    }
	  else
	    tmp = tmp->next;
	} 
      currentChrom = currentChrom->next;
      currChr++;
    }
}

long convertToBases(long totBases, int seqUnit, double value)
{
  if(seqUnit == 1)
    return ceil(totBases*value);
  else
    if(seqUnit == 2)
      return ceil(totBases*value/1000.0);
    else
      if(seqUnit == 3)
	return ceil(totBases*value/1000000.0);
      else
	return -1;
}

char JC69RBase(gsl_rng * r, char currBase)
{
  char bases[4] = {'a','c','g','t'};
  char bases_a[3] = {'c','g','t'};
  char bases_c[3] = {'a','g','t'};
  char bases_g[3] = {'a','c','t'};
  char bases_t[3] = {'a','c','g'};

  if(currBase == 'n')
    return(bases[gsl_rng_uniform_int(r, 4)]);
  else
    if(currBase == 'a')
      return(bases_a[gsl_rng_uniform_int(r, 3)]);
    else
      if(currBase == 'c')
	return(bases_c[gsl_rng_uniform_int(r, 3)]);
      else
	if(currBase == 'g')
	  return(bases_g[gsl_rng_uniform_int(r, 3)]);
	else
	  if(currBase == 't')
	    return(bases_t[gsl_rng_uniform_int(r, 3)]);
	  else
	    {
	      fprintf(stderr,"Error: base %c unknown!",currBase);
	      exit(1);
	    }
}

void getBits(unsigned int value, unsigned int noSamples, unsigned int* result)
{
  unsigned int c;
  unsigned int displaymask = 1 << (noSamples - 1);

  for(c=1; c <= noSamples; ++c)
    {
      result[c] = value & displaymask ? 1 : 0;
      value <<= 1;
    }
}

char** simulateSequences(mutation* mutation_list, int totBases, int noSamples, gsl_rng * r)
{
  unsigned int* ancPops = (unsigned int*) malloc(sizeof(unsigned int)*(noSamples+1));
  char** sequences = malloc(sizeof(char*)*noSamples);
  mutation* curr_mutList;
  long currPos=0;
  char newBase;
  
  curr_mutList = mutation_list;
  for(int i=0; i<noSamples; i++)
    sequences[i] = (char*)malloc(sizeof(char)*totBases);
  /* generate base sequences */
  for(int i=0; i<noSamples; i++)
    for(int j=0; j<totBases; j++)
      if(i==0)
	sequences[i][j] = JC69RBase(r,'n');
      else
	sequences[i][j] = sequences[0][j];
  /* add mutations */
  while(curr_mutList != NULL)
    {
      currPos = POS2BASE(totBases,curr_mutList->location);
      getBits(curr_mutList->abits,noSamples,ancPops);
      newBase = JC69RBase(r,sequences[1][currPos]);
      for(int i=1; i <= noSamples; i++)
	  if(ancPops[i])
	      sequences[i-1][currPos] = newBase;
      curr_mutList = curr_mutList->next;
    }
  return sequences;
}
