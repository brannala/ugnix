#include<unity.h>
#include<uGnix.h>

char version[] = "test";
char* test1;                                                                                                                 
char* test2;                                                                                                                 

void setUp(void) {
test1 = malloc(sizeof(char)*30);                                                                                             
test2 = malloc(sizeof(char)*30);                                                                                             
strcpy(test1,"test");                                                                                                       
strcpy(test2,"test");                                                                                                        
  // set stuff up here
}

void tearDown(void) {
    // clean stuff up here
}

void test_cstring_cmp(void) {
  TEST_ASSERT_EQUAL(0, cstring_cmp(&test1,&test2));
  
}

void test_isMissing(void) {
  TEST_ASSERT_EQUAL(1, isMissing("0"));
  TEST_ASSERT_EQUAL(1, isMissing("?")); 
  TEST_ASSERT_EQUAL(1, isMissing(".")); 
}

// not needed when using generate_test_runner.rb
int main(void) {
    UNITY_BEGIN();
    RUN_TEST(test_isMissing);
    RUN_TEST(test_cstring_cmp);
    return UNITY_END();
}
