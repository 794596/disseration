/*! \file gmp_ext.h
 */

#ifndef GMP_EXT_H
#define GMP_EXT_H

#include<gmp.h>

/**
 * Function mpz_setlimbn.
 * Sets a_n of a = (a_k ... a_0) to l.
 * If n is greater than k,
 * a is increased to size n+1 and
 * the previously unused digits
 * (a_n-1 to a_k+1) are zeroed out.
 * Leading zeroes are avoided.
 *
 * @param a the multi-precision integer.
 * @param l the single-precision integer (digit).
 * @param n the index.
 */
void mpz_setlimbn (mpz_t a, mp_limb_t l, mp_size_t n)
{
  mp_size_t i;

  if ((l == 0) && (mpz_size(a) == 0))
    return;

  /* Avoid leading zero.*/
  if ((l == 0) && (n >= mpz_size(a)))
    return;

  /* Verify that enough memory is allocated.*/
  if (n >= a->_mp_alloc)
    {
      /* Avoid to continuously allocate small pieces of memory.*/
      _mpz_realloc(a,2 * n + 1);
    }

  /* Zero out previously unused digits */
  if (n > mpz_size(a))
    {
      for(i = n; i > mpz_size(a); i--)
	a->_mp_d[i] = 0;

      /* Treat as special case to avoid that i (unsigned long) falls below zero if mpz_size(a) == 0 */
      a->_mp_d[mpz_size(a)] = 0;
    }

  a->_mp_d[n] = l;

  /* Adjust size of a.*/

  /* Eliminate leading zeroes. */
  if ((l == 0) && (mpz_size(a) > 0) && (n == mpz_size(a) - 1))
    {
      /* Check n to assure that the value of i falls not below zero. */
      if (n == 0)
	{
	  a->_mp_size = 0;
	  return;
	}

      i = n;

      while ((i-- > 0) && (a->_mp_d[i] == 0));

      if ((i == 0) && (a->_mp_d[i] == 0))
	a->_mp_size = 0;
      else
	a->_mp_size = (mpz_sgn(a) != -1) ? i+1 : -(i+1);
      return;
    }
  else
    if (mpz_size(a) < n+1)
    {
      a->_mp_size = (mpz_sgn(a) != -1) ? n+1 : -(n+1);
    }
}

/**
 * Function mp_add_limb.
 * Adds two digits and returns a two-digit result.
 *
 * @param s pointer to the digit where the least significant
 *          digit of the result should be stored.
 * @param a 1st operand.
 * @param b 2nd operand.
 * @return the carry (i.e., 2nd digit of the result).
 */
extern inline mp_limb_t mp_add_limb(mp_limb_t* s, mp_limb_t a, mp_limb_t b) { return mpn_add_n(s, &a, &b, 1); }

/**
 * Function mp_sub_limb.
 * subtracts two digits and returns a two-digit result.
 *
 * @param s pointer to the digit where the least significant
 *          digit of the result should be stored.
 * @param a 1st operand.
 * @param b 2nd operand.
 * @return the carry (i.e., 2nd digit of the result).
 */
extern inline mp_limb_t mp_sub_limb(mp_limb_t* s, mp_limb_t a, mp_limb_t b) { return mpn_sub_n(s, &a, &b, 1); }

/**
 * Function mp_mul_limb.
 * Multiplies two digits and returns a two-digit result.
 *
 * @param p pointer to the digit where the least significant
 *          digit of the result should be stored.
 * @param a 1st operand.
 * @param b 2nd operand.
 * @return the carry (i.e., 2nd digit of the result).
 */
extern inline mp_limb_t mp_mul_limb(mp_limb_t* p, mp_limb_t a, mp_limb_t b) { return mpn_mul_1(p, &a, 1, b); }

#endif /* GMP_EXT_H */
