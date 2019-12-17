#ifndef MPA_ALGS_H
#define MPA_ALGS_H

#include <gmp.h>


/**
 * Function mpz_school_add.
 * Computes the sum c of two multi-precision
 * integers a and b.
 *
 * @param c the multi-precision sum a + b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_school_add(mpz_t c, mpz_t a, mpz_t b);


void mpz_school_add_neg(mpz_t c, mpz_t a, mpz_t b);
/**
 * Function mpz_school_sub.
 * Computes the sum c of two multi-precision
 * integers a and b.
 *
 * @param c the multi-precision sum a - b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_school_sub(mpz_t c, mpz_t a, mpz_t b);

/**
 * Function mpz_mul_limb.
 * Computes the product c of a multi-precision integer a
 * and a single-precision integer b.
 *
 * @param c the multi-precision product a * b.
 * @param a the multi-precision operand.
 * @param b the single-precision operand.
 */
void mpz_mul_limb(mpz_t c, mpz_t a, mp_limb_t b);

/**
 * Function mpz_mul_base.
 * Computes the product c of a multi-precision integer a
 * and the i-th power of the base B.
 *
 * @param c the multi-precision product a * B^i.
 * @param a the multi-precision operand.
 * @param i the base exponent.
 */
void mpz_mul_base(mpz_t c, mpz_t a, mp_size_t i);

/**
 * Function mpz_school_mul.
 * Computes the product c of two multi-precision
 * integers a and b using the schoolbook method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_school_mul(mpz_t c, mpz_t a, mpz_t b);

/**
 * Function mpz_rec_mul.
 * Computes the product c of two multi-precision
 * integers a and b using a recursive multiplication method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_rec_mul(mpz_t c, mpz_t a, mpz_t b);

/**
 * Function mpz_karatsuba_mul.
 * Computes the product c of two multi-precision
 * integers a and b using Karatsuba's multiplication method.
 *
 * @param c the multi-precision product a * b.
 * @param a the first multi-precision operand.
 * @param b the second multi-precision operand.
 */
void mpz_karatsuba_mul(mpz_t c, mpz_t a, mpz_t b);

void karat_multiplying(mpz_t result, mpz_t number1, mpz_t number2);

void mpz_karatsuba_mul_par(mpz_t c, mpz_t a, mpz_t b);

void mpz_kara_mul_par(mpz_t c, mpz_t a, mpz_t b);

void multiply_kara(mpz_t c, mpz_t a, mpz_t b);

void mpz_toom3_mul(mpz_t c, mpz_t a, mpz_t b);

void mpz_toom3_mul_par(mpz_t c, mpz_t a, mpz_t b);

void mpz_toom_mul_par(mpz_t c, mpz_t a, mpz_t b);

void prnt();

#endif /* MPA_ALGS_H */
