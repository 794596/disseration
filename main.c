#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>
#include <gmp.h>
#include <string.h>
#include <math.h>
#include "multiplication.h"
#include "numArray.h"
#include "writeResults.h"

int main(int argc, char *argv[])
{
  struct timespec time1, time2;
  long long length = 10;
  double elapsed, elapsed_sec, elapsed_nano;

  //initialising variables
  mpz_t number1, number2,result,answer;
  mpz_init(result);
  mpz_init(answer);

  while (length<1000000){
    printf("%lld\n", length);
    char * num1 = file_to_array(length, "number1.txt");
    char * num2 = file_to_array(length, "number2.txt");



    // printf("%s\n", num1);
    // printf("%s\n", num2);

    mpz_init_set_str(number1, num1, 10);
    mpz_init_set_str(number2, num2, 10);
    // mpz_init_set_str(number1, num1, 10);
    // mpz_init_set_str(number2, num2, 10);
    // mpz_out_str(stdout,10,number1);
    // printf("%s\n", "");
    // mpz_out_str(stdout,10,number2);
    // printf("%s\n", "");

    // mpz_init(TrueResult);
    mpz_mul(answer, number1, number2);
    // mpz_out_str(stdout,10,answer);
    // printf("%s\n", " Answer");
    //get time = time1
    clock_gettime(CLOCK_REALTIME, &time1);


    //calling multiply to get the result of adding the two numbers
     mpz_school_mul(result, number1, number2);

     //get time = time2
     clock_gettime(CLOCK_REALTIME, &time2);

     //elapsed = (time2 - time1) * 1000
     elapsed_sec = (double)((time2.tv_sec)-(time1.tv_sec))*1000; //milliseconds
     //elapsed_nano = (time2 - time1) / 1000000
     elapsed_nano = (double)((time2.tv_nsec)-(time1.tv_nsec))/1000000;//milliseconds
     //
     if (elapsed_nano < 0) elapsed = elapsed + 1000.0;
     // elapsed = elapsed = elapsed_nano
     elapsed = elapsed_sec + elapsed_nano;
     // prints time elapsed
     if (mpz_cmp(result,answer)==0){
       writeResults(length, elapsed, "school.csv");
       printf("%s:%f\n", "school", elapsed );
       // printf("%s:%f\n", "second", elapsed_sec );
       // printf("%s:%f\n", "nano", elapsed_nano );
     } else {
       printf("%s\n", "not equal");
       writeResults(length, 0, "school.csv");
     }

     // mpz_out_str(stdout,10,result);
     // printf("%s\n", "");
     // mpz_out_str(stdout,10,answer);
     // printf("%s\n", "");


     //get time = time1
     clock_gettime(CLOCK_REALTIME, &time1);


    //calling multiply to get the result of adding the two numbers
     mpz_karatsuba_mul(result, number1, number2);




     //get time = time2
     clock_gettime(CLOCK_REALTIME, &time2);
     //elapsed = (time2 - time1) * 1000
     elapsed_sec = (double)((time2.tv_sec)-(time1.tv_sec))*1000; //milliseconds
     //elapsed_nano = (time2 - time1) / 1000000
     elapsed_nano = (double)((time2.tv_nsec)-(time1.tv_nsec))/1000000;//milliseconds
     if(elapsed_nano < 0) elapsed = elapsed + 1000.0;
     // elapsed = elapsed = elapsed_nano
     elapsed = elapsed_sec + elapsed_nano;
     // prints time elapsed
     if (mpz_cmp(result,answer)==0){
       writeResults(length, elapsed, "karatsuba.csv");
       printf("%s:%f\n", "karatsuba", elapsed );
       // printf("%s:%f\n", "second", elapsed_sec );
       // printf("%s:%f\n", "nano", elapsed_nano );
     } else {
       printf("%s\n", "not equal");
       writeResults(length, 0, "karatsuba.csv");
     }

     //get time = time1
     clock_gettime(CLOCK_REALTIME, &time1);


    //calling multiply to get the result of adding the two numbers
     mpz_kara_mul_par(result, number1, number2);



     //get time = time2
     clock_gettime(CLOCK_REALTIME, &time2);

     //elapsed = (time2 - time1) * 1000
     elapsed_sec = (double)((time2.tv_sec)-(time1.tv_sec))*1000; //milliseconds
     //elapsed_nano = (time2 - time1) / 1000000
     elapsed_nano = (double)((time2.tv_nsec)-(time1.tv_nsec))/1000000;//milliseconds
     if (elapsed_nano < 0) {
       elapsed = elapsed + 1000;
     }
     // elapsed = elapsed = elapsed_nano
     elapsed = elapsed_sec + elapsed_nano;
     // prints time elapsed
     if (mpz_cmp(result,answer)==0){
       writeResults(length, elapsed, "parallel.csv");
       printf("%s:%f\n", "parallel karatsuba", elapsed );
       printf("%s:%f\n", "second", elapsed_sec );
       printf("%s:%f\n", "nano", elapsed_nano );
     } else {
       printf("%s\n", "not equal");
       writeResults(length, 0, "parallel.csv");
     }
     //
     // mpz_out_str(stdout,10,result);
     // printf("%s\n", "");

     // printf("%d\n", (int)log10((double)length));
     length = length + 10;

   }

  return 0;

}
