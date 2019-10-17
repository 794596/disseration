#include "multiplication.h"
#include "gmp_ext.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <omp.h>


/**
 * Function mpz_school_add.
 * Computes the sum c of two multi-precision
 * integers a and b.
 *
 * @param c the multi-precision sum a + b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_school_add(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result;
    mpz_init(result);
    mp_limb_t carry = 0, carry_tmp, ai, bi, res_tmp;

    // max of len(a) and len(b)
    mp_size_t length = (mpz_size(a) > mpz_size(b)) ? mpz_size(a) : mpz_size(b);

    // Add
    for (mp_size_t i = 0; i < length; i++) {
        //get both limbs
        ai = mpz_getlimbn(a, i);
        bi = mpz_getlimbn(b, i);

        // res_tmp = ai + carry
        carry_tmp = mp_add_limb(&res_tmp, ai, carry);
        // res_tmp = res_tmp + bi
        carry = mp_add_limb(&res_tmp, res_tmp, bi);
        // carry = carry + carry_tmp
        mp_add_limb(&carry, carry, carry_tmp);
        // ci = res_tmp
        mpz_setlimbn(result, res_tmp, i);
    }
    mpz_setlimbn(result, carry, length);
    mpz_set(c, result);
    mpz_clear(result);
}


/**
 * Function mpz_mul_limb.
 * Computes the product c of a multi-precision integer a
 * and a single-precision integer b.
 *
 * @param c the multi-precision product a * b.
 * @param a the multi-precision operand.
 * @param b the single-precision operand.
 */
void mpz_mul_limb(mpz_t c, mpz_t a, mp_limb_t b) {
    mpz_t result;
    mpz_init(result);
    mp_limb_t carry = 0, carry_tmp, ai, tmp;

    for (mp_size_t i = 0; i < mpz_size(a); i++) {
        ai = mpz_getlimbn(a, i);
        carry_tmp = carry;
        //multiply
        carry = mp_mul_limb(&tmp, ai, b);
        //add the old carry
        carry_tmp = mp_add_limb(&tmp, tmp, carry_tmp);
        // (new-)carry = carry + carry_tmp
        mp_add_limb(&carry, carry, carry_tmp);
        mpz_setlimbn(result, tmp, i);
    }
    mpz_setlimbn(result, carry, mpz_size(a));
    mpz_set(c, result);
    mpz_clear(result);
}


/**
 * Function mpz_mul_base.
 * Computes the product c of a multi-precision integer a
 * and the i-th power of the base B.
 *
 * @param c the multi-precision product a * B^i.
 * @param a the multi-precision operand.
 * @param i the base exponent.
 */
void mpz_mul_base(mpz_t c, mpz_t a, mp_size_t i) {
    mpz_t result;
    mpz_init(result);

    for (mp_size_t j = 0; j < i; j++) {
        mpz_setlimbn(result, 0, j);
    }

    // copy into result
    for (mp_size_t j = 0; j < mpz_size(a); j++) {
        mpz_setlimbn(result, mpz_getlimbn(a, j), i + j);
    }
    mpz_set(c, result);
    mpz_clear(result);
}


/**
 * Function mpz_school_mul.
 * Computes the product c of two multi-precision
 * integers a and b using the schoolbook method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_school_mul(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result, tmp;
    mpz_inits(result, tmp, NULL);

    // look over all limbs of a
    for (mp_size_t i = 0; i < mpz_size(a); i++) {
        mpz_mul_base(tmp, b, i);
        mpz_mul_limb(tmp, tmp, mpz_getlimbn(a, i));
        mpz_school_add(result, result, tmp);
    }

    mpz_set(c, result);
    mpz_clears(result, tmp, NULL);
}

/**
 * Function mpz_karatsuba_mul.
 * Computes the product c of two multi-precision
 * integers a and b using Karatsuba's multiplication method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_karatsuba_mul(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result;
    mpz_init(result);


    // max of len(a) and len(b)
    mp_size_t n = (mpz_size(a) > mpz_size(b)) ? mpz_size(a) : mpz_size(b);
    if (n > 2) {
        //Prepare splitting
        // get last bit by bitwise%s and with 1
        mp_size_t m = (n & 1) ? n / 2 + 1 : n / 2;
        mpz_t a1, a0, b1, b0, a1b1, a0b0, a1a0b1b0;
        mpz_inits(a1, a0, b1, b0, a1b1, a0b0, a1a0b1b0, NULL);

        // Get upper half
        for (int i = m; i < n; ++i) {
            mpz_setlimbn(a1, mpz_getlimbn(a, i), i - m);
            mpz_setlimbn(b1, mpz_getlimbn(b, i), i - m);
        }

        // Get lower half
        for (int i = 0; i < m; ++i) {
            mpz_setlimbn(a0, mpz_getlimbn(a, i), i);
            mpz_setlimbn(b0, mpz_getlimbn(b, i), i);
        }

        // prepare multiply
        mpz_t tmp_res_1, tmp_res_2;
        mpz_inits(tmp_res_1, tmp_res_2,  NULL);



        // Now multiply and save in result
        // result = a1 * b1 * B^(2m)
        mpz_karatsuba_mul(a0b0, a0, b0);
        mpz_karatsuba_mul(a1b1, a1, b1);
        mpz_mul_base(tmp_res_1, a1b1, 2 * m);
        mpz_set(result, tmp_res_1);
        mpz_school_add(tmp_res_1, a1, a0);
        mpz_school_add(tmp_res_2, b1, b0);
        mpz_karatsuba_mul(a1a0b1b0, tmp_res_1, tmp_res_2);

        // result = result + ((a1+a0)*(b1+b0)-a1*b1-a0*b0)*B^m
        mpz_sub(tmp_res_1, a1a0b1b0, a1b1);
        mpz_sub(tmp_res_1, tmp_res_1, a0b0);
        mpz_mul_base(tmp_res_1, tmp_res_1, m);
        mpz_school_add(result, result, tmp_res_1);

        // result = result + a0*b0
        mpz_school_add(result, result,a0b0);

        // clear the variables
        mpz_clears(a1, a0, b1, b0, tmp_res_1, tmp_res_2, a1b1, a0b0, NULL);
    } else if (n > 1){
      mpz_school_mul(result,a,b);
    } else {
        // If n <= 1
        mp_limb_t carry, res;

        if (mpz_size(a) > 1 || mpz_size(b) > 1) {
            printf("Error!! Length mismatch");
        } else {
            //single precision mult
            carry = mp_mul_limb(&res, mpz_getlimbn(a, 0), mpz_getlimbn(b, 0));
            //result will consist of (carry, result)
            mpz_setlimbn(result, res, 0);
            mpz_setlimbn(result, carry, 1);
        }
    }

    // Copy result to c
    mpz_set(c, result);
    mpz_clear(result);
}

void mpz_karatsuba_mul_par(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result;
    mpz_init(result);
    // max of len(a) and len(b)
    mp_size_t n = (mpz_size(a) > mpz_size(b)) ? mpz_size(a) : mpz_size(b);
    if (n > 1) {
        //Prepare splitting
        // get last bit by bitwise and with 1
        mp_size_t m = (n & 1) ? n / 2 + 1 : n / 2;
        mpz_t a1, a0, b1, b0, a1b1, a0b0, a1a0b1b0;
        mpz_inits(a1, a0, b1, b0, a1b1, a0b0, a1a0b1b0, NULL);



        // Get upper half
        for (int i = m; i < n; ++i) {
            mpz_setlimbn(a1, mpz_getlimbn(a, i), i - m);
            mpz_setlimbn(b1, mpz_getlimbn(b, i), i - m);
        }

        // Get lower half
        for (int i = 0; i < m; ++i) {
            mpz_setlimbn(a0, mpz_getlimbn(a, i), i);
            mpz_setlimbn(b0, mpz_getlimbn(b, i), i);
        }

        // prepare multiply
        mpz_t tmp_res_1, tmp_res_2;
        mpz_inits(tmp_res_1, tmp_res_2,  NULL);

        mpz_school_add(tmp_res_1, a1, a0);
        mpz_school_add(tmp_res_2, b1, b0);

        // Now multiply and save in result
        // result = a1 * b1 * B^(2m)
          #pragma omp task shared(a0b0)
          {
            mpz_karatsuba_mul(a0b0, a0, b0);
          }
          #pragma omp task shared(a1b1)
          {
            mpz_karatsuba_mul(a1b1, a1, b1);
          }
          #pragma omp task shared(a1a0b1b0)
          {
            mpz_karatsuba_mul(a1a0b1b0, tmp_res_1, tmp_res_2);
          }
          #pragma omp taskwait
        mpz_mul_base(tmp_res_1, a1b1, 2 * m);
        mpz_set(result, tmp_res_1);


        // result = result + ((a1+a0)*(b1+b0)-a1*b1-a0*b0)*B^m
        mpz_sub(tmp_res_1, a1a0b1b0, a1b1);
        mpz_sub(tmp_res_1, tmp_res_1, a0b0);
        mpz_mul_base(tmp_res_1, tmp_res_1, m);
        mpz_school_add(result, result, tmp_res_1);

        // result = result + a0*b0
        mpz_school_add(result, result,a0b0);

        // clear the variables
        mpz_clears(a1, a0, b1, b0, tmp_res_1, tmp_res_2, a1b1, a0b0, NULL);
    } else {
        // If n <= 1
        mp_limb_t carry, res;

        if (mpz_size(a) > 1 || mpz_size(b) > 1) {
            printf("Error!! Length mismatch");
        } else {
            //single precision mult
            carry = mp_mul_limb(&res, mpz_getlimbn(a, 0), mpz_getlimbn(b, 0));
            //result will consist of (carry, result)
            mpz_setlimbn(result, res, 0);
            mpz_setlimbn(result, carry, 1);
        }
    }

    // Copy result to c
    mpz_set(c, result);
    mpz_clear(result);
}

void mpz_kara_mul_par(mpz_t c, mpz_t a, mpz_t b) {
  #pragma omp parallel num_threads(35)
  {
    #pragma omp master
    {
      mpz_karatsuba_mul_par(c,a,b);
    }
  }
}



/**
 * Function mpz_toom3_mul.
 * Computes the product c of two multi-precision
 * integers a and b using toom-3 multiplication method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_toom3_mul(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result, tmp_res_1, tmp_res_2;
    mpz_inits(result, tmp_res_1, tmp_res_2, NULL);

    // max of len(a) and len(b)
    mp_size_t n = (mpz_size(a) > mpz_size(b)) ? mpz_size(a) : mpz_size(b);
    // printf("%ld\n", mpz_size(a));
    // printf("%ld\n", mpz_size(b));
    if (n > 2) {
        //Prepare splitting
        // get last bit by bitwise and with 1
        long m;

        switch (n%3){
          case 0:
            m = n/3;
            break;
          case 1:
            m = (n+2)/3;
            break;
          case 2:
            m = (n+1)/3;
            break;
        }

        printf("m: %ld\n", m);
        mpz_t a2, a1, a0, b2, b1, b0;
        mpz_inits(a2, a1, a0, b2, b1, b0, NULL);

        // Get upper half
        for (int i = 0; i < n-m-m; i++) {
            printf("2: %d\n", i);
            mpz_setlimbn(a2, mpz_getlimbn(a, i), i);
            mpz_setlimbn(b2, mpz_getlimbn(b, i), i);
            mpz_out_str(stdout, 10, a2);
            printf("%s\n", "");
        }

        // Get lower half
        for (int i = n-m-m; i < n-m; i++) {
            printf("1: %d\n", i);
            mpz_setlimbn(a1, mpz_getlimbn(a, i), i-(n-(m+m)));
            mpz_setlimbn(b1, mpz_getlimbn(b, i), i-(n-(m+m)));
            mpz_out_str(stdout, 10, a1);
            printf("%s\n", "");
        }

        // Get lower half
        for (int i = n-m; i < n; i++) {
            printf("0: %d\n", i);
            mpz_setlimbn(a0, mpz_getlimbn(a, i), i-(n-m));
            mpz_setlimbn(b0, mpz_getlimbn(b, i), i-(n-m));
            mpz_out_str(stdout, 10, a0);
            printf("%s\n", "");
        }

        // prepare multiply
        mpz_t p, p0, p1, pm1, pm2, pinf, q, q0, q1, qm1, qm2, qinf;
        mpz_inits(p, p0, p1, pm1, pm2, pinf, q, q0, q1, qm1, qm2, qinf, NULL);

        mpz_add(p, a0,a2);
        mpz_set(p0, a0);
        mpz_add(p1, p, a1);
        mpz_sub(pm1, p, a1);
        mpz_add(tmp_res_1, pm1, a1);
        mpz_mul_ui(pm2, tmp_res_1, 2);
        mpz_sub(pm2, pm2, a0);
        mpz_set(pinf, a2);

        mpz_add(q, b0,b2);
        mpz_set(q0, b0);
        mpz_add(q1, q, b1);
        mpz_sub(qm1, q, b1);
        mpz_add(tmp_res_1, qm1, b1);
        mpz_mul_ui(qm2, tmp_res_1, 2);
        mpz_sub(qm2, qm2, b0);
        mpz_set(qinf, b2);

        mpz_t r0, r1, rm1, rm2, rinf;
        mpz_inits(r0, r1, rm1, rm2, rinf, NULL);

        mpz_mul(r0, p0, q0);
        mpz_mul(r1, p1, q1);
        mpz_mul(rm1, pm1, qm1);
        mpz_mul(rm2, pm2, qm2);
        mpz_mul(rinf, pinf, qinf);

        mpz_t res0, res1, res2, res3, res4;
        mpz_inits(res0, res1, res2, res3, res4, NULL);

        mpz_set(res0, r0);
        mpz_set(res4, rinf);
        mpz_sub(tmp_res_1, rm1, r1);
        mpz_div_ui(res3, tmp_res_1, 3);
        mpz_sub(tmp_res_1, r1, rm1);
        mpz_div_ui(res1, tmp_res_1, 2);
        mpz_sub(res2, rm1, r0);
        mpz_sub(tmp_res_1, res2, res3);
        mpz_div_ui(tmp_res_1, tmp_res_1, 2);
        mpz_mul_ui(tmp_res_2, rinf, 2);
        mpz_add(res3, tmp_res_1, tmp_res_2);
        mpz_add(res2, res2, res1);
        mpz_sub(res2, res2, res4);
        mpz_sub(res1, res1, res3);



        mpz_mul_base(res4, res4, m);
        mpz_mul_base(res3, res3, m);
        mpz_mul_base(res2, res2, m);
        mpz_mul_base(res1, res1, m);
        mpz_add(result, res0, res1);
        mpz_add(result, result, res2);
        mpz_add(result, result, res3);
        mpz_add(result, result, res4);

        // clear the variables
        mpz_clears(a1, a0, b1, b0, tmp_res_1, tmp_res_2, NULL);
    } else {
        // If n <= 1
        mp_limb_t carry, res;

        if (mpz_size(a) > 1 || mpz_size(b) > 1) {
            printf("Error!! Length mismatch");
        } else {
            //single precision mult
            carry = mp_mul_limb(&res, mpz_getlimbn(a, 0), mpz_getlimbn(b, 0));

            //result will consist of (carry, result)
            mpz_setlimbn(result, res, 0);
            mpz_setlimbn(result, carry, 1);
        }
    }

    // Copy result to c
    mpz_set(c, result);
    mpz_clear(result);
}
