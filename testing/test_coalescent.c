#include "unity.h"
#include "coalescent.h"

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
unsigned int allAnc;


void setUp(void) {
  
  allAnc = ipow(2,sampleNo);
  chromSample = create_sample(10);
  tmp2 = chromSample->chrHead;
  tmp2 = tmp2->next;

  // setup for tests test_totalAncLength and test_getRecEvent
  // chromosome 1:
  // --0--|0.2|--2--|0.42|--0--|0.64|--4--|1|
  currap = tmp2->anc;
  currap->abits = 0;
  currap->position = 0.2;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 2;
  currap->position = 0.42;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 0;
  currap->position = 0.64;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 4;
  currap->position = 1;
  currap->next = NULL;
  tmp2 = tmp2->next;
  tmp2 = tmp2->next;
  // chromosome 3:
  // --16--|0.1|--0--|0.9|--32--|1|
  currap = tmp2->anc;
  currap->abits = 16;
  currap->position = 0.1;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 0;
  currap->position = 0.9;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 32;
  currap->position = 1;
  currap->next = NULL;
  tmp2 = tmp2->next;
  tmp2 = tmp2->next;
  // chromosome 5:
  // --16--|0.1|--0--|0.9|--32--|1|
  currap = tmp2->anc;
  currap->abits = 0;
  currap->position = 0.1;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 64;
  currap->position = 0.9;
  currap->next = malloc(sizeof(ancestry));
  currap = currap->next;
  currap->abits = 0;
  currap->position = 1;
  currap->next = NULL;

  // Data for all ancestral TRUE test of TestMRCAForAll
  chrSanc = malloc(sizeof(chrsample));

  chrSanc->chrHead = malloc(sizeof(chromosome));
  chrSanc->chrHead->next = NULL;
  chrSanc->chrHead->anc = malloc(sizeof(ancestry));
  chrSanc->chrHead->anc->abits = allAnc;
  chrSanc->chrHead->anc->position = 0.5;
  chrSanc->chrHead->anc->next = malloc(sizeof(ancestry));
  chrSanc->chrHead->anc->abits = 0;
  chrSanc->chrHead->anc->position = 1.0;
  chrSanc->chrHead->anc->next = NULL;
  currC = malloc(sizeof(chromosome));
  chrSanc->chrHead->next = currC;
  currC->anc = malloc(sizeof(ancestry));
  currC->anc->abits = 0;
  currC->anc->position = 0.5;
  currC->anc->next = malloc(sizeof(ancestry));
  currC->anc->abits = allAnc;
  currC->anc->position = 1.0;
  currC->anc->next = NULL;

    // Data for all ancestral FALSE test of TestMRCAForAll
  chrSnonanc = malloc(sizeof(chrsample));
  chrSnonanc->chrHead = malloc(sizeof(chromosome));
  chrSnonanc->chrHead->next = NULL;
  chrSnonanc->chrHead->anc = malloc(sizeof(ancestry));
  chrSnonanc->chrHead->anc->abits = 8;
  chrSnonanc->chrHead->anc->position = 0.5;
  chrSnonanc->chrHead->anc->next = malloc(sizeof(ancestry));
  chrSnonanc->chrHead->anc->abits = 0;
  chrSnonanc->chrHead->anc->position = 1.0;
  chrSnonanc->chrHead->anc->next = NULL;
  currC2 = malloc(sizeof(chromosome));
  chrSnonanc->chrHead->next = currC2;
  currC2->anc = malloc(sizeof(ancestry));
  currC2->anc->abits = 0;
  currC2->anc->position = 0.5;
  currC2->anc->next = malloc(sizeof(ancestry));
  currC2->anc->abits = 9;
  currC2->anc->position = 1.0;
  currC2->anc->next = NULL;

}

void tearDown(void) {
  //  delete_sample(chromSample->chrHead);
  // free(chromSample);
  // delete_sample(chrSanc->chrHead);
  // free(chrSanc);
  // delete_sample(chrSnonanc->chrHead);
  // free(chrSnonanc);
}

void test_unionAnc(void) {
  TEST_ASSERT_EQUAL_UINT32(15, unionAnc(3,12));
  TEST_ASSERT_EQUAL_UINT32(15, unionAnc(7,15)); 
  TEST_ASSERT_EQUAL_UINT32(31, unionAnc(0,31)); 
}

void test_getChrPtr(void) {
  tmp = chromSample->chrHead;
  for(int i=0; i<3; i++)
    tmp = tmp->next;
  TEST_ASSERT_EQUAL_PTR(tmp, getChrPtr(3,chromSample));
}

void test_totalAncLength(void) {
  TEST_ASSERT_EQUAL_FLOAT(8.58, totalAncLength(chromSample));
}  

void test_getRecEvent(void) {
  getRecEvent(chromSample,2.73,&recEv);
  tmp3 = chromSample->chrHead;
  for(int i=0; i<3; i++)
    tmp3 = tmp3->next;
  TEST_ASSERT_EQUAL_PTR(tmp3, recEv.chrom);
  TEST_ASSERT_EQUAL_FLOAT(0.95,recEv.location);
}  

void test_TestMRCAForAll(void) {
  TEST_ASSERT_EQUAL(1, TestMRCAForAll(chrSanc,allAnc));
  TEST_ASSERT_EQUAL(0, TestMRCAForAll(chrSnonanc,allAnc));
}

void test_recombination(void) {
  getRecEvent(chromSample,1.42,&recEv);
  recombination(&noChrom,recEv,chromSample);
  tmp3 = chromSample->chrHead;
  for(int i=0; i<(noChrom-2); i++)
    tmp3 = tmp3->next;
  TEST_ASSERT_EQUAL_UINT32(0,tmp3->anc->abits);
  TEST_ASSERT_EQUAL_FLOAT(0.2, tmp3->anc->position);
  TEST_ASSERT_EQUAL_UINT32(2,tmp3->anc->next->abits);
  TEST_ASSERT_EQUAL_FLOAT(0.42, tmp3->anc->next->position);
  tmp3 = tmp3->next; 
  TEST_ASSERT_EQUAL_UINT32(0,tmp3->anc->abits);
  TEST_ASSERT_EQUAL_FLOAT(0.84, tmp3->anc->position);
  TEST_ASSERT_EQUAL_UINT32(4,tmp3->anc->next->abits);
  TEST_ASSERT_EQUAL_FLOAT(1.0, tmp3->anc->next->position);
  TEST_ASSERT_EQUAL_UINT32(11, noChrom);
}

void test_mergeChr(void) {
  chromosome* chr1;
  chromosome* chr2;
  chromosome* chrm;
  chr1 = malloc(sizeof(chromosome));
  chr2 = malloc(sizeof(chromosome));
  // chr1: --1--|0.2|--3--|0.7|--0--|1.0|
  chr1->anc = malloc(sizeof(ancestry));
  chr1->anc->abits = 1;
  chr1->anc->position = 0.2;
  chr1->anc->next = malloc(sizeof(ancestry));
  chr1->anc->next->abits = 3;
  chr1->anc->next->position = 0.7;
  chr1->anc->next->next = malloc(sizeof(ancestry));
  chr1->anc->next->next->abits = 0;
  chr1->anc->next->next->position = 1;
  chr1->anc->next->next->next = NULL;
  // chr2: --2--|0.5|--6--|1.0|
  chr2->anc = malloc(sizeof(ancestry));
  chr2->anc->abits = 2;
  chr2->anc->position = 0.5;
  chr2->anc->next = malloc(sizeof(ancestry));
  chr2->anc->next->abits = 6;
  chr2->anc->next->position = 1;
  // chrm: --3--|0.5|--7--|0.7|--6--|1.0|
  chrm = mergeChr(chr1,chr2);
  combineIdentAdjAncSegs(chrm);  
  TEST_ASSERT_EQUAL_UINT32(3, chrm->anc->abits);
  TEST_ASSERT_EQUAL_FLOAT(0.5, chrm->anc->position);
  TEST_ASSERT_EQUAL_UINT32(7, chrm->anc->next->abits);
  TEST_ASSERT_EQUAL_FLOAT(0.7, chrm->anc->next->position);
  TEST_ASSERT_EQUAL_UINT32(6, chrm->anc->next->next->abits);
  TEST_ASSERT_EQUAL_FLOAT(1.0, chrm->anc->next->next->position);
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
