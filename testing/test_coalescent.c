#include "unity.h"
#include "coalescent.h"

chromosome* tmp;
chrsample* chromSample;

void setUp(void) {
  chromSample = create_sample(20);
  tmp = chromSample->chrHead;
  for(int i=0; i<10; i++)
    tmp = tmp->next;

}

void tearDown(void) {
  delete_sample(chromSample->chrHead);
}

void test_unionAnc(void) {
  TEST_ASSERT_EQUAL_UINT32(15, unionAnc(3,12));
  TEST_ASSERT_EQUAL_UINT32(15, unionAnc(7,15)); 
  TEST_ASSERT_EQUAL_UINT32(31, unionAnc(0,31)); 
}

void test_getChrPtr(void) {
  TEST_ASSERT_EQUAL_PTR(tmp, getChrPtr(10,chromSample));
}

void test_totalAncLength(void) {
  TEST_ASSERT_EQUAL_FLOAT(2.0, 2.0);
}  

// not needed when using generate_test_runner.rb
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_unionAnc);
    RUN_TEST(test_getChrPtr);
    RUN_TEST(test_totalAncLength);
    return UNITY_END();
}
