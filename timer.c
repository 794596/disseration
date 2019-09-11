#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <gmp.h>
#include <string.h>
#include "multiplication.h"
#include "numArray.h"

int main(int argc, char *argv[])
{
  char * a = argv[1];

  long length = atol(a);

  struct rusage runtime;
  struct timeval t1,t2;
  //initialising variables
  mpz_t number1, number2,result;
  mpz_init(result);
  char * num1 = file_to_array(length, "number1.txt");
  char * num2 = file_to_array(length, "number2.txt");



  printf("%s\n", num1);
  printf("%s\n", num2);

  mpz_init_set_str(number1, num1, 10);
  mpz_init_set_str(number2, num2, 10);
  // mpz_init_set_str(number1, num1, 10);
  // mpz_init_set_str(number2, num2, 10);
  mpz_out_str(stdout,10,number1);
  printf("%s\n", "");
  mpz_out_str(stdout,10,number2);
  printf("%s\n", "");

  getrusage(RUSAGE_SELF, &runtime);
  t1=(runtime.ru_stime);


  //calling multiply to get the result of adding the two numbers
   mpz_karatsuba_mul(result, number1, number2);



  getrusage(RUSAGE_SELF, &runtime);
  t2=(runtime.ru_stime);

  printf("%ld\n", t1.tv_usec);
  printf("%ld\n", t2.tv_usec);

  printf("Took %ld microseconds\n", (t2.tv_usec-t1.tv_usec));
  return 0;
}
