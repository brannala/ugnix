#ifndef STDLIBS
#define STDLIBS
#include<uGnix.h>
#include<coalescent.h>
#include<unity.h>
#endif

char version[] = "test_coalescent";
chromosome* tmp=NULL;
chromosome* tmp2=NULL;
chromosome* tmp3=NULL;
chromosome* currC2;
chromosome* currC;
ancestry* currap=NULL;
chrsample* chromSample=NULL;
chrsample* chrSanc=NULL;
chrsample* chrSnonanc=NULL;
recombination_event recEv;
int sampleNo = 10;
unsigned int noChrom = 10;
bitarray* allAnc;

/* Helper to create bitarray with specific bits set */
static bitarray* make_bits(unsigned long val) {
  bitarray* ba = bitarray_create(g_noSamples);
  ba->bits[0] = val;
  ba->is_zero = (val == 0) ? 1 : 0;  /* maintain is_zero cache */
  return ba;
}

void setUp(void) {

  g_noSamples = sampleNo;
  allAnc = bitarray_full(sampleNo);
  chromSample = create_sample(10);
  tmp2 = getChrPtr(1, chromSample);

  // setup for tests test_totalAncLength and test_getRecEvent
  // chromosome 1:
  // --0--|0.2|--2--|0.42|--0--|0.64|--4--|1|
  currap = tmp2->anc;
  bitarray_clear_all(currap->abits);  // abits = 0
  currap->position = 0.2;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(2);
  currap = currap->next;
  currap->position = 0.42;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(0);
  currap = currap->next;
  currap->position = 0.64;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(4);
  currap = currap->next;
  currap->position = 1;
  currap->next = NULL;

  tmp2 = getChrPtr(3, chromSample);
  // chromosome 3:
  // --16--|0.1|--0--|0.9|--32--|1|
  currap = tmp2->anc;
  bitarray_clear_all(currap->abits);
  currap->abits->bits[0] = 16;
  currap->abits->is_zero = 0;  /* fix is_zero after direct bit manipulation */
  currap->position = 0.1;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(0);
  currap = currap->next;
  currap->position = 0.9;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(32);
  currap = currap->next;
  currap->position = 1;
  currap->next = NULL;

  tmp2 = getChrPtr(5, chromSample);
  // chromosome 5:
  // --0--|0.1|--64--|0.9|--0--|1|
  currap = tmp2->anc;
  bitarray_clear_all(currap->abits);
  currap->position = 0.1;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(64);
  currap = currap->next;
  currap->position = 0.9;
  currap->next = malloc(sizeof(ancestry));
  currap->next->abits = make_bits(0);
  currap = currap->next;
  currap->position = 1;
  currap->next = NULL;

  // Data for all ancestral TRUE test of TestMRCAForAll
  chrSanc = malloc(sizeof(chrsample));
  chrSanc->chrs = malloc(2 * sizeof(chromosome*));
  chrSanc->count = 2;
  chrSanc->capacity = 2;
  chrSanc->ancLength = 0;

  chrSanc->chrs[0] = malloc(sizeof(chromosome));
  chrSanc->chrs[0]->anc = malloc(sizeof(ancestry));
  chrSanc->chrs[0]->anc->abits = bitarray_copy(allAnc);
  chrSanc->chrs[0]->anc->position = 0.5;
  chrSanc->chrs[0]->anc->next = malloc(sizeof(ancestry));
  chrSanc->chrs[0]->anc->next->abits = make_bits(0);
  chrSanc->chrs[0]->anc->next->position = 1.0;
  chrSanc->chrs[0]->anc->next->next = NULL;

  chrSanc->chrs[1] = malloc(sizeof(chromosome));
  chrSanc->chrs[1]->anc = malloc(sizeof(ancestry));
  chrSanc->chrs[1]->anc->abits = make_bits(0);
  chrSanc->chrs[1]->anc->position = 0.5;
  chrSanc->chrs[1]->anc->next = malloc(sizeof(ancestry));
  chrSanc->chrs[1]->anc->next->abits = bitarray_copy(allAnc);
  chrSanc->chrs[1]->anc->next->position = 1.0;
  chrSanc->chrs[1]->anc->next->next = NULL;

  // Data for all ancestral FALSE test of TestMRCAForAll
  chrSnonanc = malloc(sizeof(chrsample));
  chrSnonanc->chrs = malloc(2 * sizeof(chromosome*));
  chrSnonanc->count = 2;
  chrSnonanc->capacity = 2;
  chrSnonanc->ancLength = 0;

  chrSnonanc->chrs[0] = malloc(sizeof(chromosome));
  chrSnonanc->chrs[0]->anc = malloc(sizeof(ancestry));
  chrSnonanc->chrs[0]->anc->abits = make_bits(8);
  chrSnonanc->chrs[0]->anc->position = 0.5;
  chrSnonanc->chrs[0]->anc->next = malloc(sizeof(ancestry));
  chrSnonanc->chrs[0]->anc->next->abits = make_bits(0);
  chrSnonanc->chrs[0]->anc->next->position = 1.0;
  chrSnonanc->chrs[0]->anc->next->next = NULL;

  chrSnonanc->chrs[1] = malloc(sizeof(chromosome));
  chrSnonanc->chrs[1]->anc = malloc(sizeof(ancestry));
  chrSnonanc->chrs[1]->anc->abits = make_bits(0);
  chrSnonanc->chrs[1]->anc->position = 0.5;
  chrSnonanc->chrs[1]->anc->next = malloc(sizeof(ancestry));
  chrSnonanc->chrs[1]->anc->next->abits = make_bits(9);
  chrSnonanc->chrs[1]->anc->next->position = 1.0;
  chrSnonanc->chrs[1]->anc->next->next = NULL;
}

void tearDown(void) {
  // Cleanup would go here
}

void test_unionAnc(void) {
  bitarray* a = make_bits(3);
  bitarray* b = make_bits(12);
  bitarray* result = unionAnc(a, b);
  TEST_ASSERT_EQUAL_UINT64(15, result->bits[0]);
  bitarray_free(a);
  bitarray_free(b);
  bitarray_free(result);

  a = make_bits(7);
  b = make_bits(15);
  result = unionAnc(a, b);
  TEST_ASSERT_EQUAL_UINT64(15, result->bits[0]);
  bitarray_free(a);
  bitarray_free(b);
  bitarray_free(result);

  a = make_bits(0);
  b = make_bits(31);
  result = unionAnc(a, b);
  TEST_ASSERT_EQUAL_UINT64(31, result->bits[0]);
  bitarray_free(a);
  bitarray_free(b);
  bitarray_free(result);
}

void test_getChrPtr(void) {
  tmp = getChrPtr(3, chromSample);
  TEST_ASSERT_NOT_NULL(tmp);
  // Verify it's the right chromosome by checking its ancestry
  TEST_ASSERT_EQUAL_UINT64(16, tmp->anc->abits->bits[0]);
}

void test_totalAncLength(void) {
  /* Update per-chromosome ancLen after test setup modifications */
  getChrPtr(1, chromSample)->ancLen = 0.58;  /* 0.22 + 0.36 */
  getChrPtr(3, chromSample)->ancLen = 0.2;   /* 0.1 + 0.1 */
  getChrPtr(5, chromSample)->ancLen = 0.8;   /* 0.8 */
  calcAncLength(chromSample);
  TEST_ASSERT_EQUAL_FLOAT(8.58, getAncLength(chromSample));
}

void test_getRecEvent(void) {
  /* Ensure per-chromosome ancLen is set correctly for binary search */
  getChrPtr(1, chromSample)->ancLen = 0.58;
  getChrPtr(3, chromSample)->ancLen = 0.2;
  getChrPtr(5, chromSample)->ancLen = 0.8;
  calcAncLength(chromSample);
  getRecEvent(chromSample, 2.73, &recEv);
  tmp3 = getChrPtr(3, chromSample);
  TEST_ASSERT_EQUAL_PTR(tmp3, recEv.chrom);
  TEST_ASSERT_EQUAL_FLOAT(0.95, recEv.location);
}

void test_TestMRCAForAll(void) {
  TEST_ASSERT_EQUAL(1, TestMRCAForAll(chrSanc, allAnc));
  TEST_ASSERT_EQUAL(0, TestMRCAForAll(chrSnonanc, allAnc));
}

void test_recombination(void) {
  /* Ensure per-chromosome ancLen is set correctly for binary search */
  getChrPtr(1, chromSample)->ancLen = 0.58;
  getChrPtr(3, chromSample)->ancLen = 0.2;
  getChrPtr(5, chromSample)->ancLen = 0.8;
  calcAncLength(chromSample);
  getRecEvent(chromSample, 1.42, &recEv);
  recombination(&noChrom, recEv, chromSample, NULL);

  // After recombination, noChrom should be 11
  TEST_ASSERT_EQUAL_UINT32(11, noChrom);
  TEST_ASSERT_EQUAL_INT(11, chromSample->count);

  // The recombined chromosome pieces should be at indices noChrom-2 and noChrom-1
  tmp3 = getChrPtr(noChrom - 2, chromSample);
  TEST_ASSERT_EQUAL_UINT64(0, tmp3->anc->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(0.2, tmp3->anc->position);
  TEST_ASSERT_EQUAL_UINT64(2, tmp3->anc->next->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(0.42, tmp3->anc->next->position);

  tmp3 = getChrPtr(noChrom - 1, chromSample);
  TEST_ASSERT_EQUAL_UINT64(0, tmp3->anc->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(0.84, tmp3->anc->position);
  TEST_ASSERT_EQUAL_UINT64(4, tmp3->anc->next->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(1.0, tmp3->anc->next->position);
}

void test_mergeChr(void) {
  chromosome* chr1;
  chromosome* chr2;
  chromosome* chrm;
  chr1 = malloc(sizeof(chromosome));
  chr2 = malloc(sizeof(chromosome));
  // chr1: --1--|0.2|--3--|0.7|--0--|1.0|
  chr1->anc = malloc(sizeof(ancestry));
  chr1->anc->abits = make_bits(1);
  chr1->anc->position = 0.2;
  chr1->anc->next = malloc(sizeof(ancestry));
  chr1->anc->next->abits = make_bits(3);
  chr1->anc->next->position = 0.7;
  chr1->anc->next->next = malloc(sizeof(ancestry));
  chr1->anc->next->next->abits = make_bits(0);
  chr1->anc->next->next->position = 1;
  chr1->anc->next->next->next = NULL;
  // chr2: --2--|0.5|--6--|1.0|
  chr2->anc = malloc(sizeof(ancestry));
  chr2->anc->abits = make_bits(2);
  chr2->anc->position = 0.5;
  chr2->anc->next = malloc(sizeof(ancestry));
  chr2->anc->next->abits = make_bits(6);
  chr2->anc->next->position = 1;
  chr2->anc->next->next = NULL;
  // chrm: --3--|0.5|--7--|0.7|--6--|1.0|
  chrm = mergeChr(chr1, chr2);
  combineIdentAdjAncSegs(chrm);
  TEST_ASSERT_EQUAL_UINT64(3, chrm->anc->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(0.5, chrm->anc->position);
  TEST_ASSERT_EQUAL_UINT64(7, chrm->anc->next->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(0.7, chrm->anc->next->position);
  TEST_ASSERT_EQUAL_UINT64(6, chrm->anc->next->next->abits->bits[0]);
  TEST_ASSERT_EQUAL_FLOAT(1.0, chrm->anc->next->next->position);
}

void test_coalescence(void) {
  // TODO: add coalescence test
}


// not needed when using generate_test_runner.rb
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_unionAnc);
    RUN_TEST(test_getChrPtr);
    RUN_TEST(test_totalAncLength);
    RUN_TEST(test_getRecEvent);
    RUN_TEST(test_TestMRCAForAll);
    RUN_TEST(test_recombination);
    RUN_TEST(test_mergeChr);
    return UNITY_END();
}
