#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include <gmp.h>
#include <string.h>
#include "multiplication.h"
#include "numArray.h"

/**
 * Function main
 * prints time for the the product of two mpz
 * to be produced by one of two processes.
 *
 * @param argv[1] type of multiplication
 * @param argv[2] length of mpz
 */
int main(int argc, char *argv[])
{
  char * type = argv[1];
  char * a = argv[2];
  long length = atol(a);
  struct timespec time1, time2;

  // checks to make sure that type is one of the specified strings
  if ((strcmp(type, "school") == 0) || (strcmp(type, "karat") == 0))
  {
    //initialise the mpz_t
    mpz_t number1, number2,result;
    mpz_init(result);
    mpz_init_set_str(number1, file_to_array(length, "number1.txt"), 10);
    mpz_init_set_str(number2, file_to_array(length, "number2.txt"), 10);

    // if type = "karat"
    if (strcmp(type, "karat") == 0)
    {
      //get time = time1
      clock_gettime(CLOCK_REALTIME, &time1);
      // result is product of number1 and number2 using Karatsuba multiplication
      mpz_karatsuba_mul(result, number1, number2);
      // get time = time2
      clock_gettime(CLOCK_REALTIME, &time2);
    }
    // if type = "karat"
    if (strcmp(type, "school") == 0)
    {
      //get time = time1
      clock_gettime(CLOCK_REALTIME, &time1);
      // result is product of number1 and number2 using school type multiplication
      mpz_school_mul(result, number1, number2);
      //get time = time2
      clock_gettime(CLOCK_REALTIME, &time2);
    }

    //elapsed = (time2 - time1) * 1000
    double elapsed = (double)((time2.tv_sec)-(time1.tv_sec))*1000; //milliseconds
    //elapsed_nano = (time2 - time1) / 1000000
    double elapsed_nano = (double)((time2.tv_nsec)-(time1.tv_nsec))/1000000;//milliseconds
    //
    if (elapsed_nano < 0) elapsed + 1.0;
    // elapsed = elapsed = elapsed_nano
    elapsed = elapsed + elapsed_nano;
    // prints time elapsed
    printf("%f\n", elapsed );
    return 0;
  } else
  {
    //prints correct types that can be entered
    printf("%s\n", "please enter a 'karat' or 'school'.");
  }


}
