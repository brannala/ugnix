#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#endif

/* Global sample size for bitarray allocation */
int g_noSamples = 0;

/* Global target regions for sparse simulation.
 * When set, bitarray operations are skipped for segments outside target regions. */
static const target_region_set* g_target_regions = NULL;

/* Set target regions for sparse ancestry tracking */
void set_sparse_target_regions(const target_region_set* trs) {
    g_target_regions = trs;
}

/* Check if segment [seg_start, seg_end] overlaps any target region */
static inline int segment_overlaps_targets(double seg_start, double seg_end) {
    if (!g_target_regions || !g_target_regions->active) {
        return 1;  /* No sparse mode - all segments are "targets" */
    }

    /* Quick rejection using cached bounds */
    if (seg_end <= g_target_regions->min_start ||
        seg_start >= g_target_regions->max_end) {
        return 0;
    }

    /* Check each target region */
    for (int i = 0; i < g_target_regions->n_regions; i++) {
        if (seg_start < g_target_regions->regions[i].end &&
            seg_end > g_target_regions->regions[i].start) {
            return 1;
        }
    }
    return 0;
}

/* Check if chromosome has ANY segment overlapping target regions.
 * If not, the chromosome can be pruned from the simulation. */
int chromosome_overlaps_targets(const chromosome* chr) {
    if (!g_target_regions || !g_target_regions->active) {
        return 1;  /* No sparse mode - keep all chromosomes */
    }

    double prev_pos = 0.0;
    ancestry* anc = chr->anc;

    while (anc != NULL) {
        if (segment_overlaps_targets(prev_pos, anc->position)) {
            return 1;  /* Has target-overlapping segment */
        }
        prev_pos = anc->position;
        anc = anc->next;
    }

    return 0;  /* No segments overlap targets - can be pruned */
}

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

/* O(1) check if x is a power of 2 (singleton ancestry) */
int isSingleton(unsigned int x)
{
  return x && !(x & (x - 1));
}

/* O(1) array access */
chromosome* getChrPtr(int chr, chrsample* chrom)
{
  assert(chr >= 0 && chr < chrom->count);
  return chrom->chrs[chr];
}

bitarray* unionAnc(const bitarray *anc1, const bitarray *anc2)
{
  /* Handle NULL inputs from sparse simulation */
  if (!anc1 && !anc2) {
    return NULL;  /* Both NULL - no ancestry to track */
  }
  if (!anc1) {
    return bitarray_copy(anc2);  /* Only anc2 has ancestry */
  }
  if (!anc2) {
    return bitarray_copy(anc1);  /* Only anc1 has ancestry */
  }

  /* Both have ancestry - compute union */
  bitarray *result = bitarray_create(g_noSamples);
  bitarray_union_into(result, anc1, anc2);
  return result;
}

chromosome* copy_chrom(chromosome* sourceChr)
{
  chromosome* newChr;
  ancestry* currNew;
  ancestry* currOld;
  newChr = malloc(sizeof(chromosome));
  newChr->ancLen = sourceChr->ancLen;  /* copy cached ancestral length */
  currOld = sourceChr->anc;
  newChr->anc = calloc(1, sizeof(ancestry));
  int firstAnc=1;
  currNew = newChr->anc;
  double prev_position = 0.0;  /* Track segment start for sparse check */

  while(currOld != NULL)
    {
      if(!firstAnc)
	{
	  currNew->next = calloc(1, sizeof(ancestry));
	  currNew = currNew->next;
	}
      firstAnc = 0;
      currNew->position = currOld->position;

      /* Copy MRCA cache status */
      currNew->is_mrca = currOld->is_mrca;

      /* Sparse optimization: skip bitarray copy for non-target segments.
       * But preserve zero bitarrays (non-ancestral) - don't convert to NULL
       * since NULL means "ancestral but not tracked". */
      if (segment_overlaps_targets(prev_position, currOld->position)) {
        currNew->abits = bitarray_copy(currOld->abits);
      } else if (currOld->abits == NULL || !bitarray_is_zero(currOld->abits)) {
        /* Ancestral segment outside targets - use NULL optimization */
        currNew->abits = NULL;
      } else {
        /* Non-ancestral segment (zero bitarray) - keep as zero */
        currNew->abits = bitarray_create(g_noSamples);
      }

      prev_position = currOld->position;
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
       if (tmp->abits) bitarray_free(tmp->abits);
       free(tmp);
    }
}

/* O(1) swap-and-pop deletion by index */
void delete_chrom_idx(int idx, chrsample* chrom)
{
  assert(idx >= 0 && idx < chrom->count);
  chromosome* toDelete = chrom->chrs[idx];
  delete_anc(toDelete->anc);
  free(toDelete);
  /* Swap with last element and decrement count */
  chrom->count--;
  if (idx < chrom->count) {
    chrom->chrs[idx] = chrom->chrs[chrom->count];
  }
}

void delete_sample(chrsample* chrom)
{
  for (int i = 0; i < chrom->count; i++) {
    delete_anc(chrom->chrs[i]->anc);
    free(chrom->chrs[i]);
  }
  free(chrom->chrs);
}

/* Append chromosome to array, growing if needed */
void appendChrom(chrsample* chrom, chromosome* chr)
{
  if (chrom->count >= chrom->capacity) {
    chrom->capacity *= 2;
    chrom->chrs = realloc(chrom->chrs, chrom->capacity * sizeof(chromosome*));
  }
  chrom->chrs[chrom->count++] = chr;
}

/* Calculate total ancestral length and update cache.
 * In sparse mode, NULL abits means "ancestral but not tracked" - count it. */
double calcAncLength(chrsample* chrom)
{
  ancestry* tmp_anc;
  double lastPosition = 0;
  double totLength = 0;
  for (int i = 0; i < chrom->count; i++)
    {
      tmp_anc = chrom->chrs[i]->anc;
      lastPosition = 0;
      while(tmp_anc != NULL)
	{
	  if(tmp_anc->abits == NULL || !bitarray_is_zero(tmp_anc->abits))
	      totLength += tmp_anc->position - lastPosition;
	  lastPosition = tmp_anc->position;
	  tmp_anc = tmp_anc->next;
	}
    }
  chrom->ancLength = totLength;
  return totLength;
}

/* Calculate total ACTIVE ancestral length (excluding segments at MRCA).
 * This is what should be used for coalescence and mutation rates. */
double calcActiveAncLength(chrsample* chrom, const bitarray* mrca)
{
  double totLength = 0;
  for (int i = 0; i < chrom->count; i++)
    {
      totLength += calcChromAncLengthActive(chrom->chrs[i], mrca);
    }
  return totLength;
}

/* Get cached ancestral length */
double getAncLength(const chrsample* chrom)
{
  return chrom->ancLength;
}

/* Calculate ancestral length for a single chromosome.
 * In sparse mode, NULL abits means "ancestral but not tracked in detail"
 * so those segments should contribute to length for correct rate calculations.
 * Only segments with explicit zero bitarray are non-ancestral.
 * Segments that have reached MRCA (bitarray = all samples) should NOT be
 * counted, as no more coalescence events are possible there. */
double calcChromAncLength(const chromosome* chr)
{
  double length = 0;
  double lastPos = 0;
  ancestry* tmp = chr->anc;
  while (tmp != NULL) {
    if (tmp->abits == NULL || !bitarray_is_zero(tmp->abits))
      length += tmp->position - lastPos;
    lastPos = tmp->position;
    tmp = tmp->next;
  }
  return length;
}

/* Calculate ancestral length excluding segments at MRCA.
 * This is the "active" ancestry that can still coalesce.
 * Uses cached is_mrca, bitarray_is_zero_fast, and bitarray_is_full_fast. */
double calcChromAncLengthActive(const chromosome* chr, const bitarray* mrca)
{
  double length = 0;
  double lastPos = 0;
  ancestry* tmp = chr->anc;
  while (tmp != NULL) {
    double segLen = tmp->position - lastPos;
    if (tmp->abits == NULL) {
      /* Sparse mode: count as active */
      length += segLen;
    } else if (!tmp->is_mrca) {
      /* Not cached as MRCA yet - check */
      if (!bitarray_is_zero_fast(tmp->abits)) {
        /* Fast inline check: segment can only be MRCA if all bits set */
        if (!bitarray_is_full_fast(tmp->abits)) {
          length += segLen;  /* Not full, definitely active */
        } else if (!bitarray_equal(tmp->abits, mrca)) {
          length += segLen;  /* Full but not matching (shouldn't happen) */
        } else {
          ((ancestry*)tmp)->is_mrca = 1;  /* Cache MRCA status */
        }
      }
    }
    /* If is_mrca == 1, segment is at MRCA, don't count */
    lastPos = tmp->position;
    tmp = tmp->next;
  }
  return length;
}

/* Count active segments and their total length (for progress estimation).
 * Returns segment count, stores total length in *totalLen if not NULL.
 * Uses cached is_mrca, bitarray_is_zero_fast, and bitarray_is_full_fast. */
int countActiveSegments(chrsample* chrom, const bitarray* mrca, double* totalLen)
{
  int count = 0;
  double length = 0;
  for (int i = 0; i < chrom->count; i++) {
    ancestry* tmp = chrom->chrs[i]->anc;
    double lastPos = 0;
    while (tmp != NULL) {
      double segLen = tmp->position - lastPos;
      int isActive = 0;
      if (tmp->abits == NULL) {
        isActive = 1;  /* Sparse mode */
      } else if (!tmp->is_mrca) {
        if (!bitarray_is_zero_fast(tmp->abits)) {
          /* Fast inline check: not full means definitely not MRCA */
          if (!bitarray_is_full_fast(tmp->abits)) {
            isActive = 1;
          } else if (!bitarray_equal(tmp->abits, mrca)) {
            isActive = 1;
          } else {
            ((ancestry*)tmp)->is_mrca = 1;  /* Cache MRCA status */
          }
        }
      }
      if (isActive && segLen > 1e-9) {
        count++;
        length += segLen;
      }
      lastPos = tmp->position;
      tmp = tmp->next;
    }
  }
  if (totalLen) *totalLen = length;
  return count;
}

/* Update ancLength cache after coalescence */
void updateAncLengthCoal(chrsample* chrom, double removed)
{
  chrom->ancLength -= removed;
}

/* Get active ancestral length with caching.
 * Returns cached value if valid, otherwise recalculates and caches.
 * Cache is invalidated after coalescence (when MRCA segments can change).
 * Recombination and mutation don't change total active length, so cache
 * remains valid after those events. */
double getActiveAncLength(chrsample* chrom, const bitarray* mrca)
{
  if (!chrom->activeAncValid) {
    chrom->activeAncLength = calcActiveAncLength(chrom, mrca);
    chrom->activeAncValid = 1;
  }
  return chrom->activeAncLength;
}

/* Invalidate active ancestral length cache.
 * Called after coalescence events because:
 * 1. Segments merge, reducing total ancestral length
 * 2. New MRCA segments may form
 * Recombination preserves total active length (segments split but sum is same).
 * Mutation doesn't change structure at all. */
void invalidateActiveAncLength(chrsample* chrom)
{
  chrom->activeAncValid = 0;
}

/*
eventPos is the absolute position of rec event on cumulative ancestral
chromosome material (ancLength) -> getRecEvent finds the recombinant chromosome
and relative position of recombination event on that chromosome modifies recEv

Uses binary search on precomputed cumulative sums for O(n + log n) = O(n) lookup.
The O(n) is for building cumulative array; binary search is O(log n).
For high recombination rates, this is faster than the previous O(n * segments).
*/

void getRecEvent(chrsample* chrom, double eventPos, recombination_event* recEv)
{
  getRecEventActive(chrom, eventPos, recEv, NULL);
}

/* Get recombination event location, optionally skipping MRCA segments.
 * If mrca is NULL, counts all non-zero segments (original behavior).
 * If mrca is provided, skips segments where abits == mrca. */
void getRecEventActive(chrsample* chrom, double eventPos, recombination_event* recEv,
                       const bitarray* mrca)
{
  int n = chrom->count;

  /* Build cumulative sum array using active ancestry only */
  double* cumSum = alloca((n + 1) * sizeof(double));
  cumSum[0] = 0;
  for (int i = 0; i < n; i++) {
    if (mrca)
      cumSum[i + 1] = cumSum[i] + calcChromAncLengthActive(chrom->chrs[i], mrca);
    else
      cumSum[i + 1] = cumSum[i] + chrom->chrs[i]->ancLen;
  }

  /* Binary search: find largest i where cumSum[i] <= eventPos */
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = (lo + hi + 1) / 2;  /* round up to avoid infinite loop */
    if (cumSum[mid] <= eventPos)
      lo = mid;
    else
      hi = mid - 1;
  }

  /* lo is the target chromosome index */
  double cumLen = cumSum[lo];

  /* Linear scan within the target chromosome to find exact position.
   * Skip MRCA segments if mrca is provided. */
  ancestry* tmp_anc = chrom->chrs[lo]->anc;
  double lastPosition = 0;
  double currLength = cumLen;
  double nonAncestralLength = 0;

  while (tmp_anc != NULL)
    {
      int isActive = (tmp_anc->abits == NULL) ||
                     (!bitarray_is_zero(tmp_anc->abits) &&
                      (mrca == NULL || !bitarray_equal(tmp_anc->abits, mrca)));
      if (isActive)
        currLength += tmp_anc->position - lastPosition;
      else
        nonAncestralLength += tmp_anc->position - lastPosition;
      if (currLength > eventPos)
        {
          recEv->location = eventPos + nonAncestralLength - cumLen;
          recEv->chrom = chrom->chrs[lo];
          recEv->chromIdx = lo;
          return;
        }
      lastPosition = tmp_anc->position;
      tmp_anc = tmp_anc->next;
    }
}

void recombination(unsigned int* noChrom, recombination_event recEv, chrsample* chrom)
{
  chromosome* chrtmp = recEv.chrom;
  ancestry* tmp = chrtmp->anc;
  ancestry* currAnc = NULL;
  chromosome* newLeft = malloc(sizeof(chromosome));
  chromosome* newRight = malloc(sizeof(chromosome));

  newRight->anc = calloc(1, sizeof(ancestry));
  newRight->anc->abits = bitarray_create(g_noSamples);  /* zero bitarray */
  newRight->anc->position = recEv.location;
  newRight->anc->next = NULL;
  currAnc = newRight->anc;
  while (tmp != NULL)
    {
      if (recEv.location < tmp->position)
	{
	  currAnc->next = calloc(1, sizeof(ancestry));
	  currAnc = currAnc->next;
	  currAnc->next = NULL;
	  currAnc->abits = bitarray_copy(tmp->abits);
	  currAnc->position = tmp->position;
	}
      tmp = tmp->next;
    }

  newLeft->anc = calloc(1, sizeof(ancestry));
  newLeft->anc->next = NULL;
  tmp = chrtmp->anc;
  currAnc = newLeft->anc;
  int atHead = 1;
  while (tmp->position < recEv.location)
    {
      if (!atHead)
	{
	  currAnc->next = calloc(1, sizeof(ancestry));
	  currAnc = currAnc->next;
	  currAnc->next = NULL;
	}
      currAnc->abits = bitarray_copy(tmp->abits);
      currAnc->position = tmp->position;
      tmp = tmp->next;
      atHead = 0;
    }
  if (!atHead)
    {
      currAnc->next = calloc(1, sizeof(ancestry));
      currAnc = currAnc->next;
      currAnc->next = NULL;
    }
  currAnc->abits = bitarray_copy(tmp->abits);
  currAnc->position = recEv.location;
  currAnc->next = calloc(1, sizeof(ancestry));
  currAnc = currAnc->next;
  currAnc->abits = bitarray_create(g_noSamples);  /* zero bitarray */
  currAnc->position = 1.0;
  currAnc->next = NULL;

  /* O(1) deletion using stored index */
  delete_chrom_idx(recEv.chromIdx, chrom);

  /* Check for non-empty ancestry before adding each chromosome.
     Due to floating point precision in getRecEvent, recombination location
     can occasionally land at segment boundaries, creating chromosomes
     with no ancestral material. Only add chromosomes with ancestry.
     In sparse mode, NULL abits means "ancestral but not tracked" - treat as having ancestry. */
  currAnc = newLeft->anc;
  int hasAncLeft = 0;
  while (currAnc != NULL)
    {
      if (currAnc->abits == NULL || !bitarray_is_zero(currAnc->abits)) {
        hasAncLeft = 1;
        break;
      }
      currAnc = currAnc->next;
    }

  currAnc = newRight->anc;
  int hasAncRight = 0;
  while (currAnc != NULL)
    {
      if (currAnc->abits == NULL || !bitarray_is_zero(currAnc->abits)) {
        hasAncRight = 1;
        break;
      }
      currAnc = currAnc->next;
    }

  /* Count how many chromosomes we're adding.
   * Sparse optimization: also prune chromosomes with no target-overlapping segments.
   * Track new ancLen for updating total ancLength cache. */
  int added = 0;
  double newAncLen = 0;

  if (hasAncLeft && chromosome_overlaps_targets(newLeft))
    {
      newLeft->ancLen = calcChromAncLength(newLeft);  /* cache for binary search */
      newAncLen += newLeft->ancLen;
      appendChrom(chrom, newLeft);
      added++;
    }
  else
    {
      /* Free unused chromosome (no ancestry or no target overlap) */
      delete_anc(newLeft->anc);
      free(newLeft);
    }

  if (hasAncRight && chromosome_overlaps_targets(newRight))
    {
      newRight->ancLen = calcChromAncLength(newRight);  /* cache for binary search */
      newAncLen += newRight->ancLen;
      appendChrom(chrom, newRight);
      added++;
    }
  else
    {
      /* Free unused chromosome (no ancestry or no target overlap) */
      delete_anc(newRight->anc);
      free(newRight);
    }

  /* Net change: -1 (deleted) + added */
  *noChrom = *noChrom + added - 1;
}

chromosome* mergeChr(chromosome* ptrchr1, chromosome* ptrchr2)
{
  double epsilon = 1e-10;
  chromosome* commonAnc = malloc(sizeof(chromosome));
  commonAnc->anc = calloc(1, sizeof(ancestry));
  commonAnc->anc->abits = NULL;  /* Will be set below */
  commonAnc->anc->position=0;
  commonAnc->anc->next = NULL;
  ancestry* tmp = commonAnc->anc;
  ancestry* anc1 = ptrchr1->anc;
  ancestry* anc2 = ptrchr2->anc;
  double prev_position = 0.0;  /* Track segment start for sparse check */

  while((anc1 != NULL)&&(anc2 != NULL))
    {
      /* Determine segment end position */
      double seg_end;
      if((anc1->position - anc2->position) > epsilon)
	{
	  seg_end = anc2->position;
	}
      else if((anc2->position - anc1->position) > epsilon)
	{
	  seg_end = anc1->position;
	}
      else
	{
	  seg_end = anc1->position;
	}

      /* Sparse optimization: only compute ancestry for target-overlapping segments.
       * But preserve non-ancestral status (both inputs zero) - don't use NULL
       * since NULL means "ancestral but not tracked". */
      if (tmp->abits) bitarray_free(tmp->abits);
      if (segment_overlaps_targets(prev_position, seg_end)) {
        /* Segment overlaps target - compute full ancestry */
        tmp->abits = unionAnc(anc1->abits, anc2->abits);
      } else {
        /* Segment outside targets - check if result would be ancestral */
        int anc1_is_anc = (anc1->abits == NULL || !bitarray_is_zero(anc1->abits));
        int anc2_is_anc = (anc2->abits == NULL || !bitarray_is_zero(anc2->abits));
        if (anc1_is_anc || anc2_is_anc) {
          /* At least one input is ancestral - result is ancestral */
          tmp->abits = NULL;
        } else {
          /* Both inputs non-ancestral - result is non-ancestral */
          tmp->abits = bitarray_create(g_noSamples);
        }
      }
      tmp->position = seg_end;

      /* Advance input pointers */
      if((anc1->position - anc2->position) > epsilon)
	{
	  anc2 = anc2->next;
	}
      else if((anc2->position - anc1->position) > epsilon)
	{
	  anc1 = anc1->next;
	}
      else
	{
	  anc2 = anc2->next;
	  anc1 = anc1->next;
	}

      prev_position = seg_end;

      if((anc1 != NULL)&&(anc2 != NULL))
	 {
	   tmp->next = calloc(1, sizeof(ancestry));
	   tmp->next->next = NULL;
	   tmp->next->abits = NULL;
	   tmp->next->position=0;
	   tmp = tmp->next;
	 }
    }
  return(commonAnc);
} 

/* Helper to compare ancestry bitarrays, handling NULL from sparse simulation */
static inline int ancestry_abits_equal(const bitarray* a, const bitarray* b) {
    if (a == NULL && b == NULL) return 1;  /* Both outside targets */
    if (a == NULL || b == NULL) return 0;  /* One inside, one outside */
    return bitarray_equal(a, b);
}

void combineIdentAdjAncSegs(chromosome *ptrchr)
{
  ancestry* tmp;
  ancestry* tmp_del;
  tmp = ptrchr->anc;
  while(tmp->next != NULL)
    {
      if(ancestry_abits_equal(tmp->abits, tmp->next->abits))
	{
	  tmp_del = tmp->next;
	  tmp->position = tmp->next->position;
	  tmp->next = tmp->next->next;
	  if (tmp_del->abits) bitarray_free(tmp_del->abits);
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
  chromosome* ptrchr1 = getChrPtr(pair.chr1, chrom);
  chromosome* ptrchr2 = getChrPtr(pair.chr2, chrom);

  /* Use cached ancestral lengths for incremental update */
  double len1 = ptrchr1->ancLen;
  double len2 = ptrchr2->ancLen;

  chromosome* commonAnc = mergeChr(ptrchr1, ptrchr2);
  combineIdentAdjAncSegs(commonAnc);

  /* Calculate merged length and overlap */
  double lenMerged = calcChromAncLength(commonAnc);
  commonAnc->ancLen = lenMerged;  /* cache for binary search */
  double overlap = len1 + len2 - lenMerged;

  /* Update cached ancLength incrementally */
  updateAncLengthCoal(chrom, overlap);

  /* Delete higher index first to preserve lower index validity
     (swap-and-pop changes indices) */
  int idx1 = pair.chr1;
  int idx2 = pair.chr2;
  if (idx1 > idx2) {
    delete_chrom_idx(idx1, chrom);
    delete_chrom_idx(idx2, chrom);
  } else {
    delete_chrom_idx(idx2, chrom);
    delete_chrom_idx(idx1, chrom);
  }

  /* Sparse optimization: prune chromosome if it has no target-overlapping segments */
  if (chromosome_overlaps_targets(commonAnc)) {
    /* Append merged chromosome */
    appendChrom(chrom, commonAnc);
  } else {
    /* Chromosome has no target-overlapping segments - prune it */
    delete_anc(commonAnc->anc);
    free(commonAnc);
    (*noChrom)--;  /* One fewer chromosome */
  }

  /* Invalidate active length cache - coalescence changes MRCA status */
  invalidateActiveAncLength(chrom);
}

void updateCoalescentEvents(struct coalescent_events** coalescent_list,
			    struct coalescent_events** coalescent_list_tail,
			    chrsample* chromSample, double totalTime)
{
	    struct coalescent_events* newEvent = malloc(sizeof(struct coalescent_events));

	    /* Get the last chromosome (most recently coalesced) */
	    newEvent->chr = copy_chrom(chromSample->chrs[chromSample->count - 1]);
	    combineIdentAdjAncSegs(newEvent->chr);
	    newEvent->time = totalTime;
	    newEvent->next = NULL;

	    /* O(1) append using tail pointer */
	    if(*coalescent_list == NULL)
	      {
		*coalescent_list = newEvent;
		*coalescent_list_tail = newEvent;
	      }
	    else
	      {
		(*coalescent_list_tail)->next = newEvent;
		*coalescent_list_tail = newEvent;
	      }
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

struct geneTree* getGeneTree(double lower, double upper, struct coalescent_events* coalescent_list, const bitarray* mrca)
{
  struct geneTree* geneT = NULL;
  struct geneTree* currGT = NULL;
  struct coalescent_events* localCL=coalescent_list;
  ancestry* localAnc;

  #ifdef DEBUG_GENETREE
  fprintf(stderr, "DEBUG getGeneTree: lower=%.4f upper=%.4f mrca=", lower, upper);
  bitarray_print(mrca, stderr);
  fprintf(stderr, "\n");
  struct coalescent_events* debugCL = coalescent_list;
  int eventNum = 0;
  while(debugCL != NULL) {
    fprintf(stderr, "  Event %d: time=%.2f, chr->anc->abits=", eventNum++, debugCL->time);
    bitarray_print(debugCL->chr->anc->abits, stderr);
    fprintf(stderr, "\n");
    debugCL = debugCL->next;
  }
  #endif

  while((localCL != NULL)&&((currGT == NULL)||(!bitarray_equal(currGT->abits, mrca))))
    {
      localAnc = localCL->chr->anc;
      double lastPos=0.0;
      int foundAnc=0;
      while((localAnc != NULL) && !foundAnc)
	{
	  if(!(bitarray_is_singleton(localAnc->abits) || bitarray_is_zero(localAnc->abits)))
	    {
	      if((upper <= localAnc->position)&&(lower >= lastPos))
		{
		  struct geneTree* tmpGT = geneT;
		  while((tmpGT != NULL) && !foundAnc)
		    {
		      if(bitarray_equal(tmpGT->abits, localAnc->abits))
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
		      currGT->abits = bitarray_copy(localAnc->abits);
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

void addNode(bitarray* val, double time, struct tree* lroot)
{
  while(lroot->right != NULL)
    {
      if(bitarray_equal(lroot->right->abits, val))
	{
	  lroot->right->time = time;
	  return;
	}
      else
	if(bitarray_equal(lroot->left->abits, val))
	{
	  lroot->left->time = time;
	  return;
	}
      if(bitarray_intersects(lroot->right->abits, val))
	lroot = lroot->right;
      else
	lroot = lroot->left;
    }
  lroot->right = malloc(sizeof(struct tree));
  lroot->right->right = NULL;
  lroot->right->left = NULL;
  lroot->right->abits = bitarray_copy(val);
  lroot->right->time = time;
  lroot->left = malloc(sizeof(struct tree));
  lroot->left->abits = bitarray_create(g_noSamples);
  bitarray_complement(lroot->left->abits, val);
  bitarray_intersect(lroot->left->abits, lroot->abits);
  lroot->left->left = NULL;
  lroot->left->right = NULL;
}

void splitNode(struct tree* lroot)
{
  /* Find the first set bit in lroot->abits */
  int firstBit = -1;
  for (int i = 0; i < lroot->abits->nbits; i++) {
    if (bitarray_test(lroot->abits, i)) {
      firstBit = i;
      break;
    }
  }
  if (firstBit < 0) return;

  /* Check if there are other bits set besides firstBit */
  bitarray* complement = bitarray_create(g_noSamples);
  bitarray_copy_to(complement, lroot->abits);
  bitarray_clear(complement, firstBit);

  if (!bitarray_is_zero(complement))
    {
      lroot->left = malloc(sizeof(struct tree));
      lroot->right = malloc(sizeof(struct tree));
      lroot->left->left = NULL;
      lroot->left->right = NULL;
      lroot->left->abits = bitarray_singleton(g_noSamples, firstBit);
      lroot->right->right = NULL;
      lroot->right->left = NULL;
      lroot->right->abits = complement;
    }
  else
    {
      bitarray_free(complement);
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

int binaryToChrLabel(const bitarray* ba)
{
  /* Return index of first set bit + 1 (1-indexed label) */
  for (int i = 0; i < ba->nbits; i++) {
    if (bitarray_test(ba, i))
      return i + 1;
  }
  return 0;  /* No bits set */
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
      fprintf(out, "%d:%.4f", binaryToChrLabel(node->abits), parentTime);
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
      fprintf(out, "%d;", binaryToChrLabel(root->abits));
    }
}

/* Free tree memory recursively */
void freeTree(struct tree* node)
{
  if(node == NULL) return;
  freeTree(node->left);
  freeTree(node->right);
  if(node->abits) bitarray_free(node->abits);
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
    toScreen? fprintf(stderr,"%d",binaryToChrLabel(lroot->abits)) : fprintf(tree_file,"%d",binaryToChrLabel(lroot->abits));
  return;
}

int TestMRCAForAll(chrsample* chrom, const bitarray* mrca)
{
  ancestry* tmp_anc;
  for (int i = 0; i < chrom->count; i++)
    {
      tmp_anc = chrom->chrs[i]->anc;
      while (tmp_anc != NULL)
	{
	  if (!bitarray_is_zero(tmp_anc->abits) && !bitarray_equal(tmp_anc->abits, mrca))
	    return 0; // not all zeros or all ones therefore not mrca of sample
	  tmp_anc = tmp_anc->next;
	}
    }
  return 1;
}

chrsample* create_sample(int noChrom)
{
  chrsample* chromSample = malloc(sizeof(chrsample));
  /* Initial capacity with room to grow for recombination */
  chromSample->capacity = noChrom * 4;
  chromSample->chrs = malloc(chromSample->capacity * sizeof(chromosome*));
  chromSample->count = 0;
  chromSample->ancLength = noChrom;  /* Initial: each chromosome has length 1 */
  chromSample->activeAncLength = noChrom;  /* Initially all segments are active */
  chromSample->activeAncValid = 0;  /* Will be calculated on first use */

  // create array of n sampled chromosomes
  for (int i = 0; i < noChrom; i++)
    {
      chromosome* newChrom = malloc(sizeof(chromosome));
      newChrom->anc = calloc(1, sizeof(ancestry));
      newChrom->anc->next = NULL;
      newChrom->anc->abits = bitarray_singleton(g_noSamples, i);
      newChrom->anc->position = 1.0;
      newChrom->ancLen = 1.0;  /* initial: full chromosome is ancestral */
      chromSample->chrs[chromSample->count++] = newChrom;
    }
  return chromSample;
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
 

void getMRCAs(struct mrca_list** head, chrsample* chromSample, double totalTime, const bitarray* mrca)
{
  double newlower = 0;
  double newupper = 0;
  ancestry* tmp = NULL;
  // collect information on intervals and ages of unique mrca's
  int firstInt = 1;
  /* Get the last chromosome (most recently modified) */
  chromosome* currentChrom = chromSample->chrs[chromSample->count - 1];
  tmp = currentChrom->anc;
  while ((tmp->next != NULL) || firstInt)
    {
      if (firstInt)
	{
	  if (bitarray_equal(tmp->abits, mrca))
	    {
	      newlower = 0.0;
	      newupper = tmp->position;
	      addMRCAInterval(head, newlower, newupper, totalTime);
	    }
	  firstInt = 0;
	}
      else
	{
	  if (bitarray_equal(tmp->next->abits, mrca))
	    {
	      newlower = tmp->position;
	      newupper = tmp->next->position;
	      addMRCAInterval(head, newlower, newupper, totalTime);
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

/*
getMutEvent: Find mutation location using binary search on precomputed cumulative sums.
O(n + log n) = O(n) for building array + binary search.
*/
void getMutEvent(chrsample* chrom, double eventPos, mutation* mutEv, double time)
{
  getMutEventActive(chrom, eventPos, mutEv, time, NULL);
}

/* Get mutation event location, optionally skipping MRCA segments.
 * If mrca is NULL, counts all non-zero segments (original behavior).
 * If mrca is provided, skips segments where abits == mrca. */
void getMutEventActive(chrsample* chrom, double eventPos, mutation* mutEv, double time,
                       const bitarray* mrca)
{
  int n = chrom->count;

  /* Build cumulative sum array using active ancestry only */
  double* cumSum = alloca((n + 1) * sizeof(double));
  cumSum[0] = 0;
  for (int i = 0; i < n; i++) {
    if (mrca)
      cumSum[i + 1] = cumSum[i] + calcChromAncLengthActive(chrom->chrs[i], mrca);
    else
      cumSum[i + 1] = cumSum[i] + chrom->chrs[i]->ancLen;
  }

  /* Binary search: find largest i where cumSum[i] <= eventPos */
  int lo = 0, hi = n;
  while (lo < hi) {
    int mid = (lo + hi + 1) / 2;
    if (cumSum[mid] <= eventPos)
      lo = mid;
    else
      hi = mid - 1;
  }

  double cumLen = cumSum[lo];

  /* Linear scan within target chromosome, skipping MRCA segments if provided */
  ancestry* tmp_anc = chrom->chrs[lo]->anc;
  double lastPosition = 0;
  double currLength = cumLen;
  double nonAncestralLength = 0;

  while (tmp_anc != NULL)
    {
      int isActive = (tmp_anc->abits == NULL) ||
                     (!bitarray_is_zero(tmp_anc->abits) &&
                      (mrca == NULL || !bitarray_equal(tmp_anc->abits, mrca)));
      if (isActive)
        currLength += tmp_anc->position - lastPosition;
      else
        nonAncestralLength += tmp_anc->position - lastPosition;
      if (currLength > eventPos)
        {
          mutEv->location = eventPos + nonAncestralLength - cumLen;
          mutEv->age = time;
          mutEv->abits = bitarray_copy(tmp_anc->abits);
          return;
        }
      lastPosition = tmp_anc->position;
      tmp_anc = tmp_anc->next;
    }
}

void printMutations(mutation* mutation_list, long chromTotBases, int seqUnits,
		    char* baseUnit, int noSamples, const bitarray* mrca)
{
  printf("\nMutation List\n------------------------------------\n");
  mutation* tmp1 = mutation_list;
  while(tmp1 != NULL)
    {
      if(!bitarray_equal(tmp1->abits, mrca))
	{
	  if(seqUnits)
	    printf("Pos: %ld%s Age: %.2f Chrom: ",
		   convertToBases(chromTotBases,seqUnits,tmp1->location),baseUnit,tmp1->age);
	  else
	    printf("Pos: %f Age: %.2f Chrom: ",
		   tmp1->location,tmp1->age);
	  bitarray_print(tmp1->abits, stdout);
	  printf("\n");
	}
      tmp1 = tmp1->next;
    }
}

/*
 * Write mutations as VCF format.
 * Much more memory efficient than full sequences for large sample sizes.
 */
void writeVCF(mutation* mutation_list, long chromTotBases, int noSamples,
              const bitarray* mrca, FILE* out, gsl_rng* r)
{
    /* Write VCF header */
    fprintf(out, "##fileformat=VCFv4.2\n");
    fprintf(out, "##source=coalsim\n");
    fprintf(out, "##contig=<ID=chr1,length=%ld>\n", chromTotBases);
    fprintf(out, "##INFO=<ID=AGE,Number=1,Type=Float,Description=\"Mutation age in generations\">\n");
    fprintf(out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

    /* Write sample names (haplotypes) */
    for (int i = 0; i < noSamples; i++) {
        fprintf(out, "\tsample%d", i);
    }
    fprintf(out, "\n");

    /* Bases for REF/ALT - ancestral is A, derived is random from {C,G,T} */
    char alt_bases[] = "CGT";

    /* Write each mutation as a VCF record */
    mutation* mut = mutation_list;
    while (mut != NULL) {
        /* Skip if mutation is at MRCA (all samples have it = not polymorphic) */
        if (bitarray_equal(mut->abits, mrca)) {
            mut = mut->next;
            continue;
        }

        /* Calculate base position (1-based for VCF) */
        long pos = (long)(mut->location * chromTotBases) + 1;
        if (pos < 1) pos = 1;
        if (pos > chromTotBases) pos = chromTotBases;

        /* Pick random alternate base */
        char alt = alt_bases[gsl_rng_uniform_int(r, 3)];

        /* Write VCF record */
        fprintf(out, "chr1\t%ld\t.\tA\t%c\t.\tPASS\tAGE=%.2f\tGT",
                pos, alt, mut->age);

        /* Write genotypes for each sample */
        for (int i = 0; i < noSamples; i++) {
            int has_mut = bitarray_test(mut->abits, i);
            fprintf(out, "\t%d", has_mut);
        }
        fprintf(out, "\n");

        mut = mut->next;
    }
}

void printChromosomes(chrsample* chromSample, int noSamples)
{
  ancestry* tmp = NULL;

  for (int i = 0; i < chromSample->count; i++)
    {
      printf("\nChr: %d: \n", i);
      tmp = chromSample->chrs[i]->anc;
      while (tmp != NULL)
	{
	  if ((tmp->next == NULL) || (!bitarray_equal(tmp->next->abits, tmp->abits)))
	    {
	      bitarray_print(tmp->abits, stdout);
	      printf(" %lf \n", tmp->position);
	      tmp = tmp->next;
	    }
	  else
	    tmp = tmp->next;
	}
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

/* JC69 substitution: pick uniformly from the 3 other bases
   Uses static arrays to avoid repeated stack allocation */
char JC69RBase(gsl_rng * r, char currBase)
{
  static const char bases[4] = {'a','c','g','t'};
  static const char alt[4][3] = {
    {'c','g','t'},  /* alternatives to 'a' */
    {'a','g','t'},  /* alternatives to 'c' */
    {'a','c','t'},  /* alternatives to 'g' */
    {'a','c','g'}   /* alternatives to 't' */
  };

  switch(currBase) {
    case 'n': return bases[gsl_rng_uniform_int(r, 4)];
    case 'a': return alt[0][gsl_rng_uniform_int(r, 3)];
    case 'c': return alt[1][gsl_rng_uniform_int(r, 3)];
    case 'g': return alt[2][gsl_rng_uniform_int(r, 3)];
    case 't': return alt[3][gsl_rng_uniform_int(r, 3)];
    default:
      fprintf(stderr,"Error: base %c unknown!", currBase);
      exit(1);
  }
}

/* Initialize HKY model parameters
   Precomputes normalized substitution probabilities for efficiency

   HKY model: P(i->j) proportional to:
     - kappa * pi_j  for transitions (A<->G, C<->T)
     - pi_j          for transversions (A<->C, A<->T, G<->C, G<->T)
*/
void initHKYParams(hky_params_t* params, double kappa,
                   double piA, double piC, double piG, double piT)
{
  params->kappa = kappa;
  params->pi[0] = piA;  /* A */
  params->pi[1] = piC;  /* C */
  params->pi[2] = piG;  /* G */
  params->pi[3] = piT;  /* T */

  /* Base indices: A=0, C=1, G=2, T=3
     Transitions: A<->G (0<->2), C<->T (1<->3)

     For each base, precompute normalized probabilities to other bases */

  /* From A: transition to G (kappa*piG), transversions to C,T (piC, piT) */
  double sumA = kappa * piG + piC + piT;
  params->ti_prob[0] = (kappa * piG) / sumA;  /* P(A->G) */
  params->tv_prob[0][0] = piC / sumA;          /* P(A->C) */
  params->tv_prob[0][1] = piT / sumA;          /* P(A->T) */

  /* From C: transition to T (kappa*piT), transversions to A,G (piA, piG) */
  double sumC = kappa * piT + piA + piG;
  params->ti_prob[1] = (kappa * piT) / sumC;  /* P(C->T) */
  params->tv_prob[1][0] = piA / sumC;          /* P(C->A) */
  params->tv_prob[1][1] = piG / sumC;          /* P(C->G) */

  /* From G: transition to A (kappa*piA), transversions to C,T (piC, piT) */
  double sumG = kappa * piA + piC + piT;
  params->ti_prob[2] = (kappa * piA) / sumG;  /* P(G->A) */
  params->tv_prob[2][0] = piC / sumG;          /* P(G->C) */
  params->tv_prob[2][1] = piT / sumG;          /* P(G->T) */

  /* From T: transition to C (kappa*piC), transversions to A,G (piA, piG) */
  double sumT = kappa * piC + piA + piG;
  params->ti_prob[3] = (kappa * piC) / sumT;  /* P(T->C) */
  params->tv_prob[3][0] = piA / sumT;          /* P(T->A) */
  params->tv_prob[3][1] = piG / sumT;          /* P(T->G) */
}

/* HKY substitution: sample new base according to HKY model
   Transitions (A<->G, C<->T) occur at rate kappa relative to transversions */
char HKYRBase(gsl_rng * r, char currBase, const hky_params_t* params)
{
  double u = gsl_rng_uniform(r);
  int idx;

  /* Map base to index */
  switch(currBase) {
    case 'n': {
      /* For ancestral base, sample according to equilibrium frequencies */
      double cumsum = 0;
      for (int i = 0; i < 4; i++) {
        cumsum += params->pi[i];
        if (u < cumsum) {
          static const char bases[4] = {'a','c','g','t'};
          return bases[i];
        }
      }
      return 't';  /* fallback due to floating point */
    }
    case 'a': idx = 0; break;
    case 'c': idx = 1; break;
    case 'g': idx = 2; break;
    case 't': idx = 3; break;
    default:
      fprintf(stderr,"Error: base %c unknown!", currBase);
      exit(1);
  }

  /* Sample substitution using precomputed probabilities */
  /* Order: transition first, then transversions */
  if (u < params->ti_prob[idx]) {
    /* Transition */
    static const char ti_target[4] = {'g','t','a','c'};  /* A->G, C->T, G->A, T->C */
    return ti_target[idx];
  }
  u -= params->ti_prob[idx];

  if (u < params->tv_prob[idx][0]) {
    /* First transversion */
    static const char tv1_target[4] = {'c','a','c','a'};  /* A->C, C->A, G->C, T->A */
    return tv1_target[idx];
  }

  /* Second transversion */
  static const char tv2_target[4] = {'t','g','t','g'};  /* A->T, C->G, G->T, T->G */
  return tv2_target[idx];
}

char** simulateSequences(mutation* mutation_list, int totBases, int noSamples,
                         gsl_rng * r, subst_model_t model, const hky_params_t* hky)
{
  char** sequences = malloc(sizeof(char*) * noSamples);
  mutation* curr_mutList;
  long currPos;
  char newBase;

  /* Allocate all sequences */
  for (int i = 0; i < noSamples; i++)
    sequences[i] = (char*)malloc(sizeof(char) * totBases);

  /* Generate random ancestral sequence (sample 0) */
  for (int j = 0; j < totBases; j++) {
    if (model == SUBST_HKY)
      sequences[0][j] = HKYRBase(r, 'n', hky);
    else
      sequences[0][j] = JC69RBase(r, 'n');
  }

  /* Copy ancestral to all samples using memcpy */
  for (int i = 1; i < noSamples; i++)
    memcpy(sequences[i], sequences[0], totBases);

  /* Apply mutations: use ancestral base, mutate samples with derived allele */
  curr_mutList = mutation_list;
  while (curr_mutList != NULL)
    {
      currPos = POS2BASE(totBases, curr_mutList->location);
      /* Get new base from ancestral state (sequences[0]) */
      if (model == SUBST_HKY)
        newBase = HKYRBase(r, sequences[0][currPos], hky);
      else
        newBase = JC69RBase(r, sequences[0][currPos]);
      /* Apply to all samples that carry this mutation (bit set in abits) */
      for (int i = 0; i < noSamples; i++)
        if (bitarray_test(curr_mutList->abits, i))
          sequences[i][currPos] = newBase;
      curr_mutList = curr_mutList->next;
    }
  return sequences;
}

/*
 * Target region support for sparse coalescent simulation
 */

target_region_set* create_target_regions(int n_regions, double region_length,
                                         double first_start, double spacing)
{
    if (n_regions <= 0 || region_length <= 0) return NULL;

    target_region_set* trs = malloc(sizeof(target_region_set));
    if (!trs) return NULL;

    trs->regions = malloc(n_regions * sizeof(target_region));
    if (!trs->regions) {
        free(trs);
        return NULL;
    }

    trs->n_regions = n_regions;
    trs->active = 1;

    double pos = first_start;
    double actual_spacing = spacing;

    /* If spacing is 0, distribute evenly across [0, 1] */
    if (spacing <= 0 && n_regions > 1) {
        /* Leave room for regions at both ends */
        actual_spacing = (1.0 - region_length - 2 * first_start) / (n_regions - 1);
        if (actual_spacing < region_length) {
            /* Regions would overlap, just space them by region_length */
            actual_spacing = region_length;
        }
    } else if (spacing <= 0) {
        actual_spacing = 0;  /* Single region case */
    }

    for (int i = 0; i < n_regions; i++) {
        trs->regions[i].start = pos;
        trs->regions[i].end = pos + region_length;
        /* Clamp to [0, 1] */
        if (trs->regions[i].end > 1.0) {
            trs->regions[i].end = 1.0;
        }
        pos += actual_spacing;
    }

    /* Cache bounds for fast rejection testing */
    trs->min_start = trs->regions[0].start;
    trs->max_end = trs->regions[n_regions - 1].end;

    return trs;
}

target_region_set* load_target_regions(const char* filename)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) return NULL;

    /* First pass: count lines */
    int n_regions = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        n_regions++;
    }

    if (n_regions == 0) {
        fclose(fp);
        return NULL;
    }

    target_region_set* trs = malloc(sizeof(target_region_set));
    if (!trs) {
        fclose(fp);
        return NULL;
    }

    trs->regions = malloc(n_regions * sizeof(target_region));
    if (!trs->regions) {
        free(trs);
        fclose(fp);
        return NULL;
    }

    /* Second pass: read regions */
    rewind(fp);
    int i = 0;
    while (fgets(line, sizeof(line), fp) && i < n_regions) {
        if (line[0] == '#' || line[0] == '\n') continue;
        double start, end;
        if (sscanf(line, "%lf %lf", &start, &end) == 2) {
            trs->regions[i].start = start;
            trs->regions[i].end = end;
            i++;
        }
    }

    trs->n_regions = i;
    trs->active = 1;

    /* Cache bounds for fast rejection testing */
    if (i > 0) {
        trs->min_start = trs->regions[0].start;
        trs->max_end = trs->regions[0].end;
        for (int j = 1; j < i; j++) {
            if (trs->regions[j].start < trs->min_start)
                trs->min_start = trs->regions[j].start;
            if (trs->regions[j].end > trs->max_end)
                trs->max_end = trs->regions[j].end;
        }
    }

    fclose(fp);
    return trs;
}

void free_target_regions(target_region_set* trs)
{
    if (!trs) return;
    free(trs->regions);
    free(trs);
}

/*
 * Check if a position (in 0-1 scale) falls within any target region.
 * Used to filter mutations to only target regions in sparse simulation.
 * Returns 1 if position is in a target region, 0 otherwise.
 * If regions is NULL or not active, returns 1 (accept all positions).
 */
int position_in_target_regions(double pos, const target_region_set* trs)
{
    /* If no target regions, accept all positions */
    if (!trs || !trs->active || trs->n_regions == 0) {
        return 1;
    }

    /* Quick rejection using cached bounds */
    if (pos < trs->min_start || pos > trs->max_end) {
        return 0;
    }

    /* Check each region */
    for (int i = 0; i < trs->n_regions; i++) {
        if (pos >= trs->regions[i].start && pos <= trs->regions[i].end) {
            return 1;
        }
    }

    return 0;
}

/*
 * Test if all target regions have reached MRCA.
 *
 * Optimizations:
 * 1. Single chromosome = MRCA (no need to check segments)
 * 2. Use cached min_start/max_end for quick segment rejection
 * 3. Early exit when segment starts past all targets (segments are sorted)
 * 4. Inline overlap check for single target region (common case)
 * 5. NULL abits check before bitarray_is_zero (faster)
 */
int TestMRCAForTargetRegions(chrsample* chrom, const bitarray* mrca,
                             const target_region_set* trs)
{
    /* If no target regions or not active, fall back to full check */
    if (!trs || !trs->active || trs->n_regions == 0) {
        return TestMRCAForAll(chrom, mrca);
    }

    /* Single chromosome means MRCA reached for all regions */
    if (chrom->count == 1) {
        return 1;
    }

    const double min_start = trs->min_start;
    const double max_end = trs->max_end;
    const int n_regions = trs->n_regions;
    const target_region* regions = trs->regions;

    for (int i = 0; i < chrom->count; i++) {
        ancestry* tmp_anc = chrom->chrs[i]->anc;
        double last_pos = 0.0;

        while (tmp_anc != NULL) {
            double seg_end = tmp_anc->position;

            /* Quick rejection: segment entirely before all targets */
            if (seg_end <= min_start) {
                last_pos = seg_end;
                tmp_anc = tmp_anc->next;
                continue;
            }

            /* Early exit: segment starts after all targets (segments are sorted) */
            if (last_pos >= max_end) {
                break;
            }

            /* Check overlap with target regions */
            int overlaps = 0;
            if (n_regions == 1) {
                /* Fast path for single region (common case) */
                overlaps = (last_pos < regions[0].end && seg_end > regions[0].start);
            } else {
                /* Multiple regions: linear scan (could use binary search if many) */
                for (int j = 0; j < n_regions; j++) {
                    if (last_pos < regions[j].end && seg_end > regions[j].start) {
                        overlaps = 1;
                        break;
                    }
                }
            }

            if (overlaps) {
                /* Check if this segment has reached MRCA.
                 * NULL abits in target region shouldn't happen in sparse mode,
                 * but check first as it's faster than bitarray_is_zero. */
                if (tmp_anc->abits != NULL &&
                    !bitarray_is_zero(tmp_anc->abits) &&
                    !bitarray_equal(tmp_anc->abits, mrca)) {
                    return 0;  /* Not MRCA yet for this target region */
                }
            }

            last_pos = seg_end;
            tmp_anc = tmp_anc->next;
        }
    }

    return 1;  /* All target regions have MRCA */
}
