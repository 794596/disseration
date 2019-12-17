#include "multiplication.h"
#include "gmp_ext.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

void mpz_school_add_neg(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result, a_tmp, b_tmp;
    mpz_inits(result, a_tmp, b_tmp, NULL);
    mp_limb_t carry = 0, carry_tmp, ai, bi, res_tmp;
    int a_sgn, b_sgn;
    mp_size_t length;

    a_sgn = mpz_sgn(a);
    b_sgn = mpz_sgn(b);


    // printf("%d\n", a_sgn);
    // printf("%d\n", b_sgn);
    if (a_sgn > 0 && b_sgn > 0){
        mpz_set(a_tmp, a);
        mpz_set(b_tmp, b);    // max of len(a) and len(b)

        length = (mpz_size(a_tmp) > mpz_size(b_tmp)) ? mpz_size(a_tmp) : mpz_size(b_tmp);

            // Add
        for (mp_size_t i = 0; i < length; i++) {
            //get both limbs
            ai = mpz_getlimbn(a_tmp, i);
            bi = mpz_getlimbn(b_tmp, i);

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

    } else if (a_sgn < 0 && b_sgn < 0){
      mpz_set(a_tmp, a);
      mpz_set(b_tmp, b);    // max of len(a) and len(b)

      length = (mpz_size(a_tmp) > mpz_size(b_tmp)) ? mpz_size(a_tmp) : mpz_size(b_tmp);

          // Add
      for (mp_size_t i = 0; i < length; i++) {
          //get both limbs
          ai = mpz_getlimbn(a_tmp, i);
          bi = mpz_getlimbn(b_tmp, i);

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
      mpz_neg(result, result);

    } else if (a_sgn < 0 && b_sgn > 0) {
      mpz_abs(a_tmp, a);
      mpz_abs(b_tmp, b);

      if (mpz_cmp(a_tmp,b_tmp)>0){
        mpz_sub(result, a_tmp, b_tmp);
        mpz_neg(result, result);
      } else {
        mpz_sub(result, b_tmp, a_tmp);
      }

    } else if (a_sgn > 0 && b_sgn < 0){
        mpz_abs(a_tmp, a);
        mpz_abs(b_tmp, b);
      if (mpz_cmp(a_tmp,b_tmp)>0){

        mpz_sub(result, a_tmp, b_tmp);

      } else {
        mpz_sub(result,b_tmp,  a_tmp);
        mpz_neg(result, result);
      }

    } else if (a_sgn == 0){
      mpz_set(result, b);
    } else if (b_sgn == 0){
      mpz_set(result, a);
    }

      mpz_set(c, result);
      mpz_clear(result);

}

/**
 * Function mpz_school_sub.
 * Computes the sum c of two multi-precision
 * integers a and b.
 *
 * @param c the multi-precision sum a - b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_school_sub(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result, a_tmp, b_tmp;
    mpz_inits(result, a_tmp, b_tmp, NULL);
    mp_limb_t carry = 0, carry_tmp, ai, bi, res_tmp;

    int a_sgn = mpz_sgn(a);
    int b_sgn = mpz_sgn(b);

    mpz_abs(a_tmp, a);
    mpz_abs(b_tmp, b);

    // max of len(a) and len(b)
    mp_size_t length = (mpz_size(a_tmp) > mpz_size(b_tmp)) ? mpz_size(a_tmp) : mpz_size(b_tmp);

    if (a_sgn > 0 && b_sgn < 0) {
      mpz_school_add(result, a_tmp, b_tmp);
    } else if (a_sgn < 0 && b_sgn > 0) {
      mpz_school_add(result, a_tmp, b_tmp);
      mpz_neg(result, result);
    } else if (a_sgn ==0) {
      mpz_set(result, b);
      mpz_neg(result, result);
    } else if (b_sgn ==0) {
      mpz_set(result, a);
    } else if(a_sgn < 0 && b_sgn < 0) {
      if (mpz_cmp(a_tmp, b_tmp)>0){

        mpz_sub(result, a_tmp, b_tmp);
        mpz_neg(result, result);

      } else {

        mpz_sub(result, b_tmp, a_tmp);

      }

    } else {
      if (mpz_cmp(a_tmp, b_tmp)>0){
        // Add
        for (mp_size_t i = 0; i < length; i++) {
            //get both limbs
            ai = mpz_getlimbn(a_tmp, i);
            bi = mpz_getlimbn(b_tmp, i);

            // res_tmp = ai + carry
            carry_tmp = mp_sub_limb(&res_tmp, ai, carry);
            // res_tmp = res_tmp + bi
            carry = mp_sub_limb(&res_tmp, res_tmp, bi);
            // carry = carry + carry_tmp
            mp_sub_limb(&carry, carry, carry_tmp);
            // ci = res_tmp
            mpz_setlimbn(result, res_tmp, i);
        }
        mpz_setlimbn(result, carry, length);
      } else {
        mpz_sub(result, b, a);
        mpz_neg(result, result);
      }
    }

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
    mpz_t result, a_tmp;
    mpz_inits(result, a_tmp, NULL);
    int a_sgn;

    a_sgn = mpz_sgn(a);
    mpz_abs(a_tmp, a);
    for (mp_size_t j = 0; j < i; j++) {
        mpz_setlimbn(result, 0, j);
    }

    // copy into result
    for (mp_size_t j = 0; j < mpz_size(a_tmp); j++) {
        mpz_setlimbn(result, mpz_getlimbn(a_tmp, j), i + j);
    }

    if (a_sgn < 0){
      mpz_neg(result, result);
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
    mpz_t result, tmp, a_tmp, b_tmp;
    mpz_inits(result, tmp, a_tmp, b_tmp, NULL);
    int a_sgn, b_sgn;

    mpz_abs(a_tmp, a);
    mpz_abs(b_tmp, b);

    a_sgn = mpz_sgn(a);
    b_sgn = mpz_sgn(b);

    // look over all limbs of a
    for (mp_size_t i = 0; i < mpz_size(a_tmp); i++) {
        mpz_mul_base(tmp, b_tmp, i);
        mpz_mul_limb(tmp, tmp, mpz_getlimbn(a_tmp, i));
        mpz_school_add(result, result, tmp);
    }
    if ((a_sgn >= 0 && b_sgn < 0) || (a_sgn < 0 && b_sgn >= 0)){
      mpz_neg(result, result);
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
        mpz_school_sub(tmp_res_1, a1a0b1b0, a1b1);
        mpz_school_sub(tmp_res_1, tmp_res_1, a0b0);
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
        mpz_school_sub(tmp_res_1, a1a0b1b0, a1b1);
        mpz_school_sub(tmp_res_1, tmp_res_1, a0b0);
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

void evaluation(mpz_t a2, mpz_t a1, mpz_t a0, mpz_t p,mpz_t p0,mpz_t p1,mpz_t pm1,mpz_t pm2,mpz_t pinf) {
  mpz_t tmp_res_1, tmp_res_2;
  mpz_inits(tmp_res_1, tmp_res_2, NULL);

  mpz_add(p, a0,a2);
  mpz_set(p0, a0);
  mpz_add(p1, p, a1);
  mpz_sub(pm1, p, a1);
  mpz_add(tmp_res_1, pm1, a2);
  mpz_add(pm2, tmp_res_1, tmp_res_1);
  mpz_sub(pm2, pm2, a0);
  mpz_set(pinf, a2);
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
    mpz_t result, tmp_res_1, tmp_res_2, a_tmp, b_tmp, two;
    mpz_inits(result, tmp_res_1, tmp_res_2, a_tmp, b_tmp, two, NULL);
    mpz_init_set_str(two, "2", 10);
    int a_sgn = mpz_sgn(a);
    int b_sgn = mpz_sgn(b);

    mpz_abs(a_tmp, a);
    mpz_abs(b_tmp, b);

    mp_size_t n = (mpz_size(a_tmp) > mpz_size(b_tmp)) ? mpz_size(a_tmp) : mpz_size(b_tmp);

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

        mpz_t a2, a1, a0, b2, b1, b0;
        mpz_inits(a2, a1, a0, b2, b1, b0, NULL);

        // Get upper half
        for (int i = m+m; i < n; i++) {
            mpz_setlimbn(a2, mpz_getlimbn(a_tmp, i), i-(m+m));
            mpz_setlimbn(b2, mpz_getlimbn(b_tmp, i), i-(m+m));
        }

        // Get lower half
        for (int i = m; i < m+m; i++) {
            mpz_setlimbn(a1, mpz_getlimbn(a_tmp, i), i-((m)));
            mpz_setlimbn(b1, mpz_getlimbn(b_tmp, i), i-((m)));
        }

        // Get lower half
        for (int i = 0; i < m; i++) {
            mpz_setlimbn(a0, mpz_getlimbn(a_tmp, i), i);
            mpz_setlimbn(b0, mpz_getlimbn(b_tmp, i), i);
        }

        // prepare multiply
        mpz_t p, p0, p1, pm1, pm2, pinf, q, q0, q1, qm1, qm2, qinf;
        mpz_inits(p, p0, p1, pm1, pm2, pinf, q, q0, q1, qm1, qm2, qinf, NULL);

        evaluation(a2, a1, a0, p, p0, p1, pm1, pm2, pinf);
        evaluation(b2, b1, b0, q, q0, q1, qm1, qm2, qinf);

        mpz_t r0, r1, rm1, rm2, rinf;
        mpz_inits(r0, r1, rm1, rm2, rinf, NULL);

        mpz_toom3_mul(r0, p0, q0);
        mpz_toom3_mul(r1, p1, q1);
        mpz_toom3_mul(rm1, pm1, qm1);
        mpz_toom3_mul(tmp_res_1, pm1, qm1);
        mpz_toom3_mul(rm2, pm2, qm2);
        mpz_toom3_mul(rinf, pinf, qinf);


        mpz_t res0, res1, res2, res3, res4;
        mpz_inits(res0, res1, res2, res3, res4, NULL);

        mpz_set(res0, r0);
        mpz_set(res4, rinf);
        mpz_sub(tmp_res_1, rm2, r1);
        mpz_divexact_ui(res3, tmp_res_1, 3);
        mpz_sub(tmp_res_1, r1, rm1);
        mpz_div_2exp(res1, tmp_res_1, 1);
        mpz_sub(res2, rm1, r0);
        mpz_sub(tmp_res_1, res2, res3);
        mpz_sub(tmp_res_2, res2, res3);
        mpz_sub(tmp_res_2, tmp_res_2, tmp_res_1);

        mpz_div_2exp(tmp_res_1, tmp_res_1, 1);
        mpz_add(tmp_res_2, rinf, rinf);
        mpz_add(res3, tmp_res_1, tmp_res_2);
        mpz_add(res2, res2, res1);
        mpz_sub(res2, res2, res4);
        mpz_sub(res1, res1, res3);

        mpz_mul_base(res4, res4,m*4);
        mpz_mul_base(res3, res3,m*3);
        mpz_mul_base(res2, res2,m*2);
        mpz_mul_base(res1, res1,m*1);


        mpz_add(result, res0, res1);
        mpz_add(result, result, res2);
        mpz_add(result, result, res3);
        mpz_add(result, result, res4);

        // clear the variables
        mpz_clears(a1, a0, b1, b0, tmp_res_1, tmp_res_2, NULL);

      } else {
        // If n <= 1
        //mp_limb_t carry, res;

        if (mpz_size(a_tmp) > 2 || mpz_size(b_tmp) > 2) {

            printf("Error!! Length mismatch");

        } else {

          mpz_mul(result, a_tmp, b_tmp);
            // //single precision mult
            //carry = mp_mul_limb(&res, mpz_getlimbn(a_tmp, 0), mpz_getlimbn(b, 0));
            //
            // //result will consist of (carry, result)
            // mpz_setlimbn(result, res, 0);
            // mpz_setlimbn(result, carry, 1);
        }
    }
    if(((b_sgn > 0) && (a_sgn <= 0))||((a_sgn > 0) && (b_sgn <= 0))){

      mpz_neg(result, result);

    }

    // Copy result to c
    mpz_set(c, result);
    mpz_clear(result);
}

void mpz_toom_mul_par(mpz_t c, mpz_t a, mpz_t b) {
  mpz_t result, tmp_res_1, tmp_res_2, a_tmp, b_tmp, two;
  mpz_inits(result, tmp_res_1, tmp_res_2, a_tmp, b_tmp, two, NULL);
  mpz_init_set_str(two, "2", 10);
  int a_sgn = mpz_sgn(a);
  int b_sgn = mpz_sgn(b);

  // mpz_out_str(stdout, 10, b);
  // printf("%s\n", "");
  //
  // mpz_out_str(stdout, 10, b_tmp);
  // printf("%s\n", "");

  mpz_abs(a_tmp, a);
  mpz_abs(b_tmp, b);

  // mpz_out_str(stdout, 10, b);
  // printf("%s\n", "");
  //
  // mpz_out_str(stdout, 10, b_tmp);
  // printf("%s\n", "");

  // max of len(a) and len(b)
  mp_size_t n = (mpz_size(a_tmp) > mpz_size(b_tmp)) ? mpz_size(a_tmp) : mpz_size(b_tmp);

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

      mpz_t a2, a1, a0, b2, b1, b0;
      mpz_inits(a2, a1, a0, b2, b1, b0, NULL);




      // Get upper half
      for (int i = m+m; i < n; i++) {
          mpz_setlimbn(a2, mpz_getlimbn(a_tmp, i), i-(m+m));
          mpz_setlimbn(b2, mpz_getlimbn(b_tmp, i), i-(m+m));
      }


      // Get lower half
      for (int i = m; i < m+m; i++) {
          mpz_setlimbn(a1, mpz_getlimbn(a_tmp, i), i-((m)));
          mpz_setlimbn(b1, mpz_getlimbn(b_tmp, i), i-((m)));

      }


      // Get lower half
      for (int i = 0; i < m; i++) {
          mpz_setlimbn(a0, mpz_getlimbn(a_tmp, i), i);
          mpz_setlimbn(b0, mpz_getlimbn(b_tmp, i), i);
      }

      // prepare multiply
      mpz_t p, p0, p1, pm1, pm2, pinf, q, q0, q1, qm1, qm2, qinf;
      mpz_inits(p, p0, p1, pm1, pm2, pinf, q, q0, q1, qm1, qm2, qinf, NULL);

      evaluation(a2, a1, a0, p, p0, p1, pm1, pm2, pinf);
      evaluation(b2, b1, b0, q, q0, q1, qm1, qm2, qinf);

      mpz_t r0, r1, rm1, rm2, rinf;
      mpz_inits(r0, r1, rm1, rm2, rinf, NULL);

      #pragma omp task shared(r0)
      {
        mpz_toom_mul_par(r0, p0, q0);
      }
      #pragma omp task shared(r1)
      {
        mpz_toom_mul_par(r1, p1, q1);
      }
      #pragma omp task shared(rm1)
      {
        mpz_toom_mul_par(rm1, pm1, qm1);
      }
      #pragma omp task shared(rm2)
      {
        mpz_toom_mul_par(rm2, pm2, qm2);
      }
      #pragma omp task shared(rinf)
      {
        mpz_toom_mul_par(rinf, pinf, qinf);
      }
      #pragma omp taskwait

      mpz_t res0, res1, res2, res3, res4;
      mpz_inits(res0, res1, res2, res3, res4, NULL);

      mpz_set(res0, r0);
      mpz_set(res4, rinf);
      mpz_sub(tmp_res_1, rm2, r1);
      mpz_divexact_ui(res3, tmp_res_1, 3);
      mpz_sub(tmp_res_1, r1, rm1);
      mpz_div_2exp(res1, tmp_res_1, 1);
      mpz_sub(res2, rm1, r0);
      mpz_sub(tmp_res_1, res2, res3);
      mpz_sub(tmp_res_2, res2, res3);
      mpz_sub(tmp_res_2, tmp_res_2, tmp_res_1);

      mpz_div_2exp(tmp_res_1, tmp_res_1, 1);
      mpz_add(tmp_res_2, rinf, rinf);
      mpz_add(res3, tmp_res_1, tmp_res_2);
      mpz_add(res2, res2, res1);
      mpz_sub(res2, res2, res4);
      mpz_sub(res1, res1, res3);

      mpz_mul_base(res4, res4,m*4);
      mpz_mul_base(res3, res3,m*3);
      mpz_mul_base(res2, res2,m*2);
      mpz_mul_base(res1, res1,m*1);


      mpz_add(result, res0, res1);
      mpz_add(result, result, res2);
      mpz_add(result, result, res3);
      mpz_add(result, result, res4);

      // clear the variables
      mpz_clears(a1, a0, b1, b0, tmp_res_1, tmp_res_2, NULL);

    } else {
      // If n <= 1
      //mp_limb_t carry, res;

      if (mpz_size(a_tmp) > 2 || mpz_size(b_tmp) > 2) {

          printf("Error!! Length mismatch");

      } else {

        mpz_mul(result, a_tmp, b_tmp);
          // //single precision mult
          //carry = mp_mul_limb(&res, mpz_getlimbn(a_tmp, 0), mpz_getlimbn(b, 0));
          //
          // //result will consist of (carry, result)
          // mpz_setlimbn(result, res, 0);
          // mpz_setlimbn(result, carry, 1);
      }
  }
  if(((b_sgn > 0) && (a_sgn <= 0))||((a_sgn > 0) && (b_sgn <= 0))){

    mpz_neg(result, result);

  }

  // Copy result to c
  mpz_set(c, result);
  mpz_clear(result);
}

void mpz_toom3_mul_par(mpz_t c, mpz_t a, mpz_t b) {
  #pragma omp parallel num_threads(35)
  {
    #pragma omp master
    {
      mpz_toom_mul_par(c,a,b);
    }
  }
}
