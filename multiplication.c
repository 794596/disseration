#include "multiplication.h"
#include "gmp_ext.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>


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
        // res_tmp = res_temp + bi
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
 * Function mpz_rec_mul.
 * Computes the product c of two multi-precision
 * integers a and b using a recursive multiplication method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_rec_mul(mpz_t c, mpz_t a, mpz_t b) {
    mpz_t result;
    mpz_init(result);

    // max of len(a) and len(b)
    mp_size_t n = (mpz_size(a) > mpz_size(b)) ? mpz_size(a) : mpz_size(b);

    if (n > 1) {
        //Prepare splitting
        // get last bit by bitwise and with 1
        mp_size_t m = (n & 1) ? n / 2 + 1 : n / 2;
        mpz_t a1, a0, b1, b0;
        mpz_inits(a1, a0, b1, b0, NULL);

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
        mpz_inits(tmp_res_1, tmp_res_2, NULL);

        // Now multiply and save in result
        // result = a1 * b1 * B^(2m)
        mpz_rec_mul(tmp_res_1, a1, b1);
        mpz_mul_base(tmp_res_1, tmp_res_1, 2 * m);
        mpz_set(result, tmp_res_1);

        // result = result + (a1 * b0 + a0 * b1) * B^m
        mpz_rec_mul(tmp_res_2, a1, b0);
        mpz_rec_mul(tmp_res_1, a0, b1);
        mpz_school_add(tmp_res_2, tmp_res_2, tmp_res_1);
        mpz_mul_base(tmp_res_2, tmp_res_2, m);
        mpz_school_add(result, result, tmp_res_2);

        // result = result + a0 * b0
        mpz_rec_mul(tmp_res_1, a0, b0);
        mpz_school_add(result, result, tmp_res_1);

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

    if (n > 1) {
        //Prepare splitting
        // get last bit by bitwise and with 1
        mp_size_t m = (n & 1) ? n / 2 + 1 : n / 2;
        mpz_t a1, a0, b1, b0;
        mpz_inits(a1, a0, b1, b0, NULL);

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
        mpz_inits(tmp_res_1, tmp_res_2, NULL);

        // Now multiply and save in result
        // result = a1 * b1 * B^(2m)
        mpz_karatsuba_mul(tmp_res_1, a1, b1);
        mpz_mul_base(tmp_res_1, tmp_res_1, 2 * m);
        mpz_set(result, tmp_res_1);

        // result = result + ((a1+a0)*(b1+b0)-a1*b1-a0*b0)*B^m
        mpz_school_add(tmp_res_1, a1, a0);
        mpz_school_add(tmp_res_2, b1, b0);
        mpz_karatsuba_mul(tmp_res_1, tmp_res_1, tmp_res_2);
        mpz_karatsuba_mul(tmp_res_2, a1, b1);
        mpz_sub(tmp_res_1, tmp_res_1, tmp_res_2);
        mpz_karatsuba_mul(tmp_res_2, a0, b0);
        mpz_sub(tmp_res_1, tmp_res_1, tmp_res_2);
        mpz_mul_base(tmp_res_1, tmp_res_1, m);
        mpz_school_add(result, result, tmp_res_1);

        // result = result + a0*b0
        mpz_karatsuba_mul(tmp_res_1, a0, b0);
        mpz_school_add(result, result, tmp_res_1);

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


/**
 * Function mpz_montgomery_trans.
 * Performs montgomery transformation. Maps the operands
 * to another representation system.
 *
 * @param x_tilde the transformed multi-precision result.
 * @param x the multi-precision operand.
 */

void mpz_montgomery_trans(mpz_t x_tilde, mpz_t x, mpz_t m) {
    mpz_t result, r;
    mpz_init(result);

    mpz_init_set_ui(r, 1);
    mpz_mul_base(r, r, mpz_size(m));
    mpz_karatsuba_mul(result, x, r);
    mpz_mod(result, result, m);

    mpz_set(x_tilde, result);
    mpz_clears(result, r, NULL);
}


/**
 * Function mpz_montgomery_get_m_prime.
 * Computes the negative inverse of the modulus. The result
 * m' is needed for montgomery reduction
 *
 * @param m_prime the negative inverse of m.
 * @param m the modulus.
 */
void mpz_montgomery_get_m_prime(mpz_t m_prime, mpz_t m) {
    mpz_t invers, r;
    mpz_init(invers);
    mpz_init_set_ui(r, 1);

    mpz_mul_base(r, r, mpz_size(m));
    mpz_invert(invers, m, r);
    mpz_sub(m_prime, r, invers);

    mpz_clears(invers, r, NULL);
}


/**
 * Function mpz_montgomery_red.
 * Computes montgomery reduction of a multi-precision product
 * T in montgomery respresentation system
 *
 * @param c the reduced multi-precision product
 * @param T the input multi-precision product.
 * @param m the modulus.
 * @param m_prime the negative inverse of m.
 */
void mpz_montgomery_red(mpz_t c, mpz_t T, mpz_t m, mpz_t m_prime) {
    mpz_t result, r, u;
    mpz_inits(r, result, u, NULL);
    mpz_set_ui(r, 1);

    mpz_mul_base(r, r, mpz_size(m));
    mpz_karatsuba_mul(u, T, m_prime);

    // result = u mod r
    for (long i = 0; i < mpz_size(r) - 1; ++i) {
        mpz_setlimbn(result, mpz_getlimbn(u, i), i);
    }

    // result = result * m + T
    mpz_karatsuba_mul(result, result, m);
    mpz_school_add(result, result, T);

    for (int i = 0; i < mpz_size(result) - mpz_size(m); i++) {
        mpz_setlimbn(c, mpz_getlimbn(result, i + mpz_size(m)), i);
    }

    // Check c < m
    if (mpz_cmp(c, m) > 0) {
        mpz_sub(c, c, m);
    }

    mpz_clears(result, r, u, NULL);
}
