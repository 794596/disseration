#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <gmp.h>
#include "multiplication.h"


// long long *toom(long long m, long long n)
// {
//   long long m2, m1, m0, n2, n1, n0, length, b;
//
//   length = (log10(m) > log10(n)) ? log10(m)+1 : log10(n)+1;
//
//   switch (length%3) {
//     case 0:
//       b = (length/3);
//       break;
//     case 1:
//       b = ((length+2)/3);
//       break;
//     case 2:
//       b = ((length+1)/3);
//       break;
//   }
//
//   m0 = m % (long long)pow(10, b);
//   m1 = (m % (long long)pow(10, b*2)-m0)/(long long)pow(10, b);
//   m2 = (m-(m1*(long long)pow(10, b))-m0)/(long long)pow(10, b*2);
//
//   n0 = n % (long long)pow(10, b);
//   n1 = ((n % (long long)pow(10, b*2))-n0)/(long long)pow(10, b);
//   n2 = (n-(n1*(long long)pow(10, b))-n0)/(long long)pow(10, b*2);
//
//   long long p, p0, p1, pn1, pn2, pinf;
//
//   p = m0+m2;
//   p0 = m0;
//   p1 = p+m1;
//   pn1 = p-m1;
//   pn2 = ((pn1+m2)*2)-m0;
//   pinf = m2;
//
//   long long q, q0, q1, qn1, qn2, qinf;
//
//   q = n0+n2;
//   q0 = n0;
//   q1 = q+n1;
//   qn1 = q-n1;
//   qn2 = ((qn1+n2)*2)-n0;
//   qinf = n2;
//
//   long long  r0, r1, rn1, rn2, rinf;
//
//   r0 = p0*q0;
//   r1 = p1*q1;
//   rn1 = pn1*qn1;
//   rn2 = pn2*qn2;
//   rinf = pinf*qinf;
//
//   long long  res0, res1, res2, res3, res4, result;
//
//   res0 = r0;
//   res4 = rinf;
//   res3 = (rn2-r1)/3;
//   res1 = (r1-rn1)/2;
//   res2 = (rn1-r0);
//   res3 = ((res2-res3)/2)+(2*rinf);
//   res2 = res2+res1-res4;
//   res1 = res1-res3;
//
//   result = (res0+(res1*(long long)pow(10, b))+(res2*(long long)pow(10, b*2))+(res3*(long long)pow(10, b*3))+(res4*(long long)pow(10, b*4)));
//   return (long long *)result;
// }



int main(){
  // unsigned long long *c;
  // c = (unsigned long long *) malloc(sizeof(unsigned long long)*30);
  // unsigned long long m = 12345678912345;
  // unsigned long long n = 9876543211245;
  // printf("%lld\n", m*n);
  // c = toom(m, n);
  // printf("%lld\n", (long long)c);
  mpz_t result, result_karat, result_toom, a, b, x, y, m, n;
  mpz_inits(result, result_karat, result_toom, x, y, m, n, NULL);
  mpz_init_set_str(a, "123456789012345678901234567890123456789012345", 10);
  mpz_init_set_str(b, "987654321098765432109876543210987654321098765", 10);
  mpz_toom3_mul(result_toom, a,b);
  printf("%s", "Toom:      ");
  mpz_out_str(stdout, 16, result_toom);
  printf("\n%s", "Karatsuba: ");
  mpz_karatsuba_mul(result_karat, a, b);
  mpz_out_str(stdout, 16, result_karat);
  printf("%s\n            ", "");
  mpz_sub(result, result_karat, result_toom);
  mpz_out_str(stdout, 16, result);
  printf("%s\n", "");

  mpz_init_set_str(x, "1000000", 10);
  mpz_init_set_str(y, "100000", 10);

  mpz_school_add(result, y, x);


  mpz_out_str(stdout, 10, result);
  printf("%s\n", "");

  mpz_init_set_str(n, "-100", 10);
  mpz_init_set_str(m, "12345", 10);

  mpz_school_add(result, n, m);

  mpz_out_str(stdout, 10, result);
  printf("%s\n", "");

  mpz_set_str(n, "-10000000000000000000000000000000000000000", 10);
  mpz_set_str(m, "-12345678900000000000000000000000000000000", 10);

  mpz_school_add(result, n, m);

  mpz_out_str(stdout, 10, result);
  printf("%s\n", "");

    return 0;

}
