//
//  main.cpp
//  test
//
//  Created by ZhaoYoccio on 16/9/28.
//  Copyright © 2016年 ZhaoYoccio. All rights reserved.
//
#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

typedef long mp_size_t;
typedef unsigned long mp_bitcnt_t;
typedef unsigned long mp_limb_t;
typedef mp_limb_t *mp_ptr;
typedef const mp_limb_t *mp_srcptr;


typedef struct
{
    int _mp_alloc;		/* Number of *limbs* allocated and pointed
                         to by the _mp_d field.  */
    int _mp_size;			/* abs(_mp_size) is the number of limbs the
                             last field points to.  If _mp_size is
                             negative this is a negative number.  */
    mp_limb_t *_mp_d;		/* Pointer to the limbs.  */
} __mpz_struct;
typedef const __mpz_struct *mpz_srcptr;
typedef __mpz_struct mpz_t[1];
#define mpn_invert_limb(x) mpn_invert_3by2 ((x), 0)

//////////////////
//LINE 55
#define GMP_LIMB_BITS (sizeof(mp_limb_t) * CHAR_BIT)
//LINE 57
#define GMP_LIMB_MAX (~ (mp_limb_t) 0)
//LINE 58
#define GMP_LIMB_HIGHBIT ((mp_limb_t) 1 << (GMP_LIMB_BITS - 1))
//LINE 60
#define GMP_HLIMB_BIT ((mp_limb_t) 1 << (GMP_LIMB_BITS / 2))
//LINE 61
#define GMP_LLIMB_MASK (GMP_HLIMB_BIT - 1)
//LINE 66
#define GMP_ABS(x) ((x) >= 0 ? (x) : -(x))
//LINE 67
#define GMP_NEG_CAST(T,x) (-((T)((x) + 1) - 1))
//LINE 70
#define GMP_MAX(a, b) ((a) > (b) ? (a) : (b))
//LINE 72
#define gmp_assert_nocarry(x) do { \
mp_limb_t __cy = x;		   \
assert (__cy == 0);		   \
} while (0)
//LINE 77
#define gmp_clz(count, x) do {						\
mp_limb_t __clz_x = (x);						\
unsigned __clz_c;							\
for (__clz_c = 0;							\
(__clz_x & ((mp_limb_t) 0xff << (GMP_LIMB_BITS - 8))) == 0;	\
__clz_c += 8)							\
__clz_x <<= 8;							\
for (; (__clz_x & GMP_LIMB_HIGHBIT) == 0; __clz_c++)		\
__clz_x <<= 1;							\
(count) = __clz_c;							\
} while (0)
//LINE 96
#define gmp_add_ssaaaa(sh, sl, ah, al, bh, bl) \
do {									\
mp_limb_t __x;							\
__x = (al) + (bl);							\
(sh) = (ah) + (bh) + (__x < (al));					\
(sl) = __x;								\
} while(0)
//LINE 104
#define gmp_sub_ddmmss(sh, sl, ah, al, bh, bl) \
do {									\
mp_limb_t __x;							\
__x = (al) - (bl);							\
(sh) = (ah) - (bh) - ((al) < (bl));					\
(sl) = __x;								\
} while (0)
//LINE 112
#define gmp_umul_ppmm(w1, w0, u, v)					\
do {									\
mp_limb_t __x0, __x1, __x2, __x3;					\
unsigned __ul, __vl, __uh, __vh;					\
mp_limb_t __u = (u), __v = (v);					\
\
__ul = __u & GMP_LLIMB_MASK;					\
__uh = __u >> (GMP_LIMB_BITS / 2);					\
__vl = __v & GMP_LLIMB_MASK;					\
__vh = __v >> (GMP_LIMB_BITS / 2);					\
\
__x0 = (mp_limb_t) __ul * __vl;					\
__x1 = (mp_limb_t) __ul * __vh;					\
__x2 = (mp_limb_t) __uh * __vl;					\
__x3 = (mp_limb_t) __uh * __vh;					\
\
__x1 += __x0 >> (GMP_LIMB_BITS / 2);/* this can't give carry */	\
__x1 += __x2;		/* but this indeed can */		\
if (__x1 < __x2)		/* did we get it? */			\
__x3 += GMP_HLIMB_BIT;	/* yes, add it in the proper pos. */	\
\
(w1) = __x3 + (__x1 >> (GMP_LIMB_BITS / 2));			\
(w0) = (__x1 << (GMP_LIMB_BITS / 2)) + (__x0 & GMP_LLIMB_MASK);	\
} while (0)
//LINE 137
#define gmp_udiv_qrnnd_preinv(q, r, nh, nl, d, di)			\
do {									\
mp_limb_t _qh, _ql, _r, _mask;					\
gmp_umul_ppmm (_qh, _ql, (nh), (di));				\
gmp_add_ssaaaa (_qh, _ql, _qh, _ql, (nh) + 1, (nl));		\
_r = (nl) - _qh * (d);						\
_mask = -(mp_limb_t) (_r > _ql); /* both > and >= are OK */		\
_qh += _mask;							\
_r += _mask & (d);							\
if (_r >= (d))							\
{									\
_r -= (d);							\
_qh++;								\
}									\
\
(r) = _r;								\
(q) = _qh;								\
} while (0)
//LINE 156
#define gmp_udiv_qr_3by2(q, r1, r0, n2, n1, n0, d1, d0, dinv)		\
do {									\
mp_limb_t _q0, _t1, _t0, _mask;					\
gmp_umul_ppmm ((q), _q0, (n2), (dinv));				\
gmp_add_ssaaaa ((q), _q0, (q), _q0, (n2), (n1));			\
\
/* Compute the two most significant limbs of n - q'd */		\
(r1) = (n1) - (d1) * (q);						\
gmp_sub_ddmmss ((r1), (r0), (r1), (n0), (d1), (d0));		\
gmp_umul_ppmm (_t1, _t0, (d0), (q));				\
gmp_sub_ddmmss ((r1), (r0), (r1), (r0), _t1, _t0);			\
(q)++;								\
\
/* Conditionally adjust q and the remainders */			\
_mask = - (mp_limb_t) ((r1) >= _q0);				\
(q) += _mask;							\
gmp_add_ssaaaa ((r1), (r0), (r1), (r0), _mask & (d1), _mask & (d0)); \
if ((r1) >= (d1))							\
{									\
if ((r1) > (d1) || (r0) >= (d0))				\
{								\
(q)++;							\
gmp_sub_ddmmss ((r1), (r0), (r1), (r0), (d1), (d0));	\
}								\
}									\
} while (0)
//LINE 190
#define MP_SIZE_T_SWAP(x, y)						\
do {									\
mp_size_t __mp_size_t_swap__tmp = (x);				\
(x) = (y);								\
(y) = __mp_size_t_swap__tmp;					\
} while (0)
//LINE 202
#define MP_PTR_SWAP(x, y)						\
do {									\
mp_ptr __mp_ptr_swap__tmp = (x);					\
(x) = (y);								\
(y) = __mp_ptr_swap__tmp;						\
} while (0)
//LINE 232
#define MPZ_SRCPTR_SWAP(x, y)						\
do {									\
mpz_srcptr __mpz_srcptr_swap__tmp = (x);				\
(x) = (y);								\
(y) = __mpz_srcptr_swap__tmp;					\
} while (0)
//LINE 244
static void
gmp_die (const char *msg)
{
    fprintf (stderr, "%s\n", msg);
    abort();
}
//LINE 250
static void *
gmp_default_alloc (size_t size)
{
    void *p;
    
    assert (size > 0);
    
    p = malloc (size);
    if (!p)
        gmp_die("gmp_default_alloc: Virtual memory exhausted.");
    
    return p;
}
//LINE 264
static void *
gmp_default_realloc (void *old, size_t old_size, size_t new_size)
{
    void * p;
    
    p = realloc (old, new_size);
    
    if (!p)
        gmp_die("gmp_default_realloc: Virtual memory exhausted.");
    
    return p;
}
//LINE 277
static void
gmp_default_free (void *p, size_t size)
{
    free (p);
}
//LINE 283
static void * (*gmp_allocate_func) (size_t) = gmp_default_alloc;
//LINE 284
static void * (*gmp_reallocate_func) (void *, size_t, size_t) = gmp_default_realloc;
//LINE 285
static void (*gmp_free_func) (void *, size_t) = gmp_default_free;
//LINE 319
#define gmp_xalloc(size) ((*gmp_allocate_func)((size)))
//LINE 320
#define gmp_free(p) ((*gmp_free_func) ((p), 0))
//LINE 322
static mp_ptr
gmp_xalloc_limbs (mp_size_t size)
{
    return (mp_ptr) gmp_xalloc (size * sizeof (mp_limb_t));
}
//LINE 328
static mp_ptr
gmp_xrealloc_limbs (mp_ptr old, mp_size_t size)
{
    assert (size > 0);
    return (mp_ptr) (*gmp_reallocate_func) (old, 0, size * sizeof (mp_limb_t));
}
//LINE 338
void
mpn_copyi (mp_ptr d, mp_srcptr s, mp_size_t n)
{
    mp_size_t i;
    for (i = 0; i < n; i++)
        d[i] = s[i];
}
//LINE 353
int
mpn_cmp (mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
    while (--n >= 0)
    {
        if (ap[n] != bp[n])
            return ap[n] > bp[n] ? 1 : -1;
    }
    return 0;
}
//LINE 364
static int
mpn_cmp4 (mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    if (an != bn)
        return an < bn ? -1 : 1;
    else
        return mpn_cmp (ap, bp, an);
}
//LINE 373
static mp_size_t
mpn_normalized_size (mp_srcptr xp, mp_size_t n)
{
    while (n > 0 && xp[n-1] == 0)
        --n;
    return n;
}
//LINE 394
mp_limb_t
mpn_add_1 (mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b)
{
    mp_size_t i;
    
    assert (n > 0);
    i = 0;
    do
    {
        mp_limb_t r = ap[i] + b;
        /* Carry out */
        b = (r < b);
        rp[i] = r;
    }
    while (++i < n);
    
    return b;
}
//LINE 413
mp_limb_t
mpn_add_n (mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
    mp_size_t i;
    mp_limb_t cy;
    
    for (i = 0, cy = 0; i < n; i++)
    {
        mp_limb_t a, b, r;
        a = ap[i]; b = bp[i];
        r = a + cy;
        cy = (r < cy);
        r += b;
        cy += (r < b);
        rp[i] = r;
    }
    return cy;
}
//LINE 432
mp_limb_t
mpn_add (mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    mp_limb_t cy;
    
    assert (an >= bn);
    
    cy = mpn_add_n (rp, ap, bp, bn);
    if (an > bn)
        cy = mpn_add_1 (rp + bn, ap + bn, an - bn, cy);
    return cy;
}
//LINE 445
mp_limb_t
mpn_sub_1 (mp_ptr rp, mp_srcptr ap, mp_size_t n, mp_limb_t b)
{
    mp_size_t i;
    
    assert (n > 0);
    
    i = 0;
    do
    {
        mp_limb_t a = ap[i];
        /* Carry out */
        mp_limb_t cy = a < b;;
        rp[i] = a - b;
        b = cy;
    }
    while (++i < n);
    
    return b;
}
//LINE 466
mp_limb_t
mpn_sub_n (mp_ptr rp, mp_srcptr ap, mp_srcptr bp, mp_size_t n)
{
    mp_size_t i;
    mp_limb_t cy;
    
    for (i = 0, cy = 0; i < n; i++)
    {
        mp_limb_t a, b;
        a = ap[i]; b = bp[i];
        b += cy;
        cy = (b < cy);
        cy += (a < b);
        rp[i] = a - b;
    }
    return cy;
}
//LINE 484
mp_limb_t
mpn_sub (mp_ptr rp, mp_srcptr ap, mp_size_t an, mp_srcptr bp, mp_size_t bn)
{
    mp_limb_t cy;
    
    assert (an >= bn);
    
    cy = mpn_sub_n (rp, ap, bp, bn);
    if (an > bn)
        cy = mpn_sub_1 (rp + bn, ap + bn, an - bn, cy);
    return cy;
}
//LINE 497
mp_limb_t
mpn_mul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
    mp_limb_t ul, cl, hpl, lpl;
    
    assert (n >= 1);
    
    cl = 0;
    do
    {
        ul = *up++;
        gmp_umul_ppmm (hpl, lpl, ul, vl);
        
        lpl += cl;
        cl = (lpl < cl) + hpl;
        
        *rp++ = lpl;
    }
    while (--n != 0);
    
    return cl;
}
//LINE 520
mp_limb_t
mpn_addmul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
    mp_limb_t ul, cl, hpl, lpl, rl;
    
    assert (n >= 1);
    
    cl = 0;
    do
    {
        ul = *up++;
        gmp_umul_ppmm (hpl, lpl, ul, vl);
        
        lpl += cl;
        cl = (lpl < cl) + hpl;
        
        rl = *rp;
        lpl = rl + lpl;
        cl += lpl < rl;
        *rp++ = lpl;
    }
    while (--n != 0);
    
    return cl;
}
//LINE 546
mp_limb_t
mpn_submul_1 (mp_ptr rp, mp_srcptr up, mp_size_t n, mp_limb_t vl)
{
    mp_limb_t ul, cl, hpl, lpl, rl;
    
    assert (n >= 1);
    
    cl = 0;
    do
    {
        ul = *up++;
        gmp_umul_ppmm (hpl, lpl, ul, vl);
        
        lpl += cl;
        cl = (lpl < cl) + hpl;
        
        rl = *rp;
        lpl = rl - lpl;
        cl += lpl > rl;
        *rp++ = lpl;
    }
    while (--n != 0);
    
    return cl;
}
//LINE 572
mp_limb_t
mpn_mul (mp_ptr rp, mp_srcptr up, mp_size_t un, mp_srcptr vp, mp_size_t vn)
{
    assert (un >= vn);
    assert (vn >= 1);
    
    /* We first multiply by the low order limb. This result can be
     stored, not added, to rp. We also avoid a loop for zeroing this
     way. */
    
    rp[un] = mpn_mul_1 (rp, up, un, vp[0]);
    
    /* Now accumulate the product of up[] and the next higher limb from
     vp[]. */
    
    while (--vn >= 1)
    {
        rp += 1, vp += 1;
        rp[un] = mpn_addmul_1 (rp, up, un, vp[0]);
    }
    return rp[un];
}

//LINE 607
mp_limb_t
mpn_lshift (mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
    mp_limb_t high_limb, low_limb;
    unsigned int tnc;
    mp_limb_t retval;
    
    assert (n >= 1);
    assert (cnt >= 1);
    assert (cnt < GMP_LIMB_BITS);
    
    up += n;
    rp += n;
    
    tnc = GMP_LIMB_BITS - cnt;
    low_limb = *--up;
    retval = low_limb >> tnc;
    high_limb = (low_limb << cnt);
    
    while (--n != 0)
    {
        low_limb = *--up;
        *--rp = high_limb | (low_limb >> tnc);
        high_limb = (low_limb << cnt);
    }
    *--rp = high_limb;
    
    return retval;
}
//LINE 662
mp_limb_t
mpn_rshift (mp_ptr rp, mp_srcptr up, mp_size_t n, unsigned int cnt)
{
    mp_limb_t high_limb, low_limb;
    unsigned int tnc;
    mp_limb_t retval;
    
    assert (n >= 1);
    assert (cnt >= 1);
    assert (cnt < GMP_LIMB_BITS);
    
    tnc = GMP_LIMB_BITS - cnt;
    high_limb = *up++;
    retval = (high_limb << tnc);
    low_limb = high_limb >> cnt;
    
    while (--n != 0)
    {
        high_limb = *up++;
        *rp++ = low_limb | (high_limb << tnc);
        low_limb = high_limb >> cnt;
    }
    *rp = low_limb;
    
    return retval;
}
//LINE 706
mp_limb_t
mpn_invert_3by2 (mp_limb_t u1, mp_limb_t u0)
{
    mp_limb_t r, p, m;
    unsigned ul, uh;
    unsigned ql, qh;
    
    /* First, do a 2/1 inverse. */
    /* The inverse m is defined as floor( (B^2 - 1 - u1)/u1 ), so that 0 <
     * B^2 - (B + m) u1 <= u1 */
    assert (u1 >= GMP_LIMB_HIGHBIT);
    
    ul = u1 & GMP_LLIMB_MASK;
    uh = u1 >> (GMP_LIMB_BITS / 2);
    
    qh = ~u1 / uh;
    r = ((~u1 - (mp_limb_t) qh * uh) << (GMP_LIMB_BITS / 2)) | GMP_LLIMB_MASK;
    
    p = (mp_limb_t) qh * ul;
    /* Adjustment steps taken from udiv_qrnnd_c */
    if (r < p)
    {
        qh--;
        r += u1;
        if (r >= u1) /* i.e. we didn't get carry when adding to r */
            if (r < p)
            {
                qh--;
                r += u1;
            }
    }
    r -= p;
    
    /* Do a 3/2 division (with half limb size) */
    p = (r >> (GMP_LIMB_BITS / 2)) * qh + r;
    ql = (p >> (GMP_LIMB_BITS / 2)) + 1;
    
    /* By the 3/2 method, we don't need the high half limb. */
    r = (r << (GMP_LIMB_BITS / 2)) + GMP_LLIMB_MASK - ql * u1;
    
    if (r >= (p << (GMP_LIMB_BITS / 2)))
    {
        ql--;
        r += u1;
    }
    m = ((mp_limb_t) qh << (GMP_LIMB_BITS / 2)) + ql;
    if (r >= u1)
    {
        m++;
        r -= u1;
    }
    
    if (u0 > 0)
    {
        mp_limb_t th, tl;
        r = ~r;
        r += u0;
        if (r < u0)
        {
            m--;
            if (r >= u1)
            {
                m--;
                r -= u1;
            }
            r -= u1;
        }
        gmp_umul_ppmm (th, tl, u0, m);
        r += th;
        if (r < th)
        {
            m--;
            m -= ((r > u1) | ((r == u1) & (tl > u0)));
        }
    }
    
    return m;
}
//LINE 785
struct gmp_div_inverse
{
    /* Normalization shift count. */
    unsigned shift;
    /* Normalized divisor (d0 unused for mpn_div_qr_1) */
    mp_limb_t d1, d0;
    /* Inverse, for 2/1 or 3/2. */
    mp_limb_t di;
};
//LINE 795
static void
mpn_div_qr_1_invert (struct gmp_div_inverse *inv, mp_limb_t d)
{
    unsigned shift;
    
    assert (d > 0);
    gmp_clz (shift, d);
    inv->shift = shift;
    inv->d1 = d << shift;
    inv->di = mpn_invert_limb (inv->d1);
}
//LINE 807
static void
mpn_div_qr_2_invert (struct gmp_div_inverse *inv,
                     mp_limb_t d1, mp_limb_t d0)
{
    unsigned shift;
    
    assert (d1 > 0);
    gmp_clz (shift, d1);
    inv->shift = shift;
    if (shift > 0)
    {
        d1 = (d1 << shift) | (d0 >> (GMP_LIMB_BITS - shift));
        d0 <<= shift;
    }
    inv->d1 = d1;
    inv->d0 = d0;
    inv->di = mpn_invert_3by2 (d1, d0);
}

//LINE 826
static void
mpn_div_qr_invert (struct gmp_div_inverse *inv,
                   mp_srcptr dp, mp_size_t dn)
{
    assert (dn > 0);
    
    if (dn == 1)
        mpn_div_qr_1_invert (inv, dp[0]);
    else if (dn == 2)
        mpn_div_qr_2_invert (inv, dp[1], dp[0]);
    else
    {
        unsigned shift;
        mp_limb_t d1, d0;
        
        d1 = dp[dn-1];
        d0 = dp[dn-2];
        assert (d1 > 0);
        gmp_clz (shift, d1);
        inv->shift = shift;
        if (shift > 0)
        {
            d1 = (d1 << shift) | (d0 >> (GMP_LIMB_BITS - shift));
            d0 = (d0 << shift) | (dp[dn-3] >> (GMP_LIMB_BITS - shift));
        }
        inv->d1 = d1;
        inv->d0 = d0;
        inv->di = mpn_invert_3by2 (d1, d0);
    }
}
//LINE 859
static mp_limb_t
mpn_div_qr_1_preinv (mp_ptr qp, mp_srcptr np, mp_size_t nn,
                     const struct gmp_div_inverse *inv)
{
    mp_limb_t d, di;
    mp_limb_t r;
    mp_ptr tp = NULL;
    
    if (inv->shift > 0)
    {
        tp = gmp_xalloc_limbs (nn);
        r = mpn_lshift (tp, np, nn, inv->shift);
        np = tp;
    }
    else
        r = 0;
    
    d = inv->d1;
    di = inv->di;
    while (--nn >= 0)
    {
        mp_limb_t q;
        
        gmp_udiv_qrnnd_preinv (q, r, r, np[nn], d, di);
        if (qp)
            qp[nn] = q;
    }
    if (inv->shift > 0)
        gmp_free (tp);
    
    return r >> inv->shift;
}
//LINE 922
static void
mpn_div_qr_2_preinv (mp_ptr qp, mp_ptr rp, mp_srcptr np, mp_size_t nn,
                     const struct gmp_div_inverse *inv)
{
    unsigned shift;
    mp_size_t i;
    mp_limb_t d1, d0, di, r1, r0;
    mp_ptr tp;
    
    assert (nn >= 2);
    shift = inv->shift;
    d1 = inv->d1;
    d0 = inv->d0;
    di = inv->di;
    
    if (shift > 0)
    {
        tp = gmp_xalloc_limbs (nn);
        r1 = mpn_lshift (tp, np, nn, shift);
        np = tp;
    }
    else
        r1 = 0;
    
    r0 = np[nn - 1];
    
    i = nn - 2;
    do
    {
        mp_limb_t n0, q;
        n0 = np[i];
        gmp_udiv_qr_3by2 (q, r1, r0, r1, r0, n0, d1, d0, di);
        
        if (qp)
            qp[i] = q;
    }
    while (--i >= 0);
    
    if (shift > 0)
    {
        assert ((r0 << (GMP_LIMB_BITS - shift)) == 0);
        r0 = (r0 >> shift) | (r1 << (GMP_LIMB_BITS - shift));
        r1 >>= shift;
        
        gmp_free (tp);
    }
    
    rp[1] = r1;
    rp[0] = r0;
}
//LINE 986
static void
mpn_div_qr_pi1 (mp_ptr qp,
                mp_ptr np, mp_size_t nn, mp_limb_t n1,
                mp_srcptr dp, mp_size_t dn,
                mp_limb_t dinv)
{
    mp_size_t i;
    
    mp_limb_t d1, d0;
    mp_limb_t cy, cy1;
    mp_limb_t q;
    
    assert (dn > 2);
    assert (nn >= dn);
    
    d1 = dp[dn - 1];
    d0 = dp[dn - 2];
    
    assert ((d1 & GMP_LIMB_HIGHBIT) != 0);
    /* Iteration variable is the index of the q limb.
     *
     * We divide <n1, np[dn-1+i], np[dn-2+i], np[dn-3+i],..., np[i]>
     * by            <d1,          d0,        dp[dn-3],  ..., dp[0] >
     */
    
    i = nn - dn;
    do
    {
        mp_limb_t n0 = np[dn-1+i];
        
        if (n1 == d1 && n0 == d0)
        {
            q = GMP_LIMB_MAX;
            mpn_submul_1 (np+i, dp, dn, q);
            n1 = np[dn-1+i];	/* update n1, last loop's value will now be invalid */
        }
        else
        {
            gmp_udiv_qr_3by2 (q, n1, n0, n1, n0, np[dn-2+i], d1, d0, dinv);
            
            cy = mpn_submul_1 (np + i, dp, dn-2, q);
            
            cy1 = n0 < cy;
            n0 = n0 - cy;
            cy = n1 < cy1;
            n1 = n1 - cy1;
            np[dn-2+i] = n0;
            
            if (cy != 0)
            {
                n1 += d1 + mpn_add_n (np + i, np + i, dp, dn - 1);
                q--;
            }
        }
        
        if (qp)
            qp[i] = q;
    }
    while (--i >= 0);
    
    np[dn - 1] = n1;
}
//LINE 1049
static void
mpn_div_qr_preinv (mp_ptr qp, mp_ptr np, mp_size_t nn,
                   mp_srcptr dp, mp_size_t dn,
                   const struct gmp_div_inverse *inv)
{
    assert (dn > 0);
    assert (nn >= dn);
    
    if (dn == 1)
        np[0] = mpn_div_qr_1_preinv (qp, np, nn, inv);
    else if (dn == 2)
        mpn_div_qr_2_preinv (qp, np, np, nn, inv);
    else
    {
        mp_limb_t nh;
        unsigned shift;
        
        assert (inv->d1 == dp[dn-1]);
        assert (inv->d0 == dp[dn-2]);
        assert ((inv->d1 & GMP_LIMB_HIGHBIT) != 0);
        
        shift = inv->shift;
        if (shift > 0)
            nh = mpn_lshift (np, np, nn, shift);
        else
            nh = 0;
        
        mpn_div_qr_pi1 (qp, np, nn, nh, dp, dn, inv->di);
        
        if (shift > 0)
            gmp_assert_nocarry (mpn_rshift (np, np, dn, shift));
    }
}

//LINE 1083
static void
mpn_div_qr (mp_ptr qp, mp_ptr np, mp_size_t nn, mp_srcptr dp, mp_size_t dn)
{
    struct gmp_div_inverse inv;
    mp_ptr tp = NULL;
    
    assert (dn > 0);
    assert (nn >= dn);
    
    mpn_div_qr_invert (&inv, dp, dn);
    if (dn > 2 && inv.shift > 0)
    {
        tp = gmp_xalloc_limbs (dn);
        gmp_assert_nocarry (mpn_lshift (tp, dp, dn, inv.shift));
        dp = tp;
    }
    mpn_div_qr_preinv (qp, np, nn, dp, dn, &inv);
    if (tp)
        gmp_free (tp);
}

//LINE 1105
/* MPN base conversion. */
static unsigned
mpn_base_power_of_two_p (unsigned b)
{
    switch (b)
    {
        case 2: return 1;
        case 4: return 2;
        case 8: return 3;
        case 16: return 4;
        case 32: return 5;
        case 64: return 6;
        case 128: return 7;
        case 256: return 8;
        default: return 0;
    }
}
//LINE 1123
struct mpn_base_info
{
    /* bb is the largest power of the base which fits in one limb, and
     exp is the corresponding exponent. */
    unsigned exp;
    mp_limb_t bb;
};
//LINE 1131
static void
mpn_get_base_info (struct mpn_base_info *info, mp_limb_t b)
{
    mp_limb_t m;
    mp_limb_t p;
    unsigned exp;
    
    m = GMP_LIMB_MAX / b;
    for (exp = 1, p = b; p <= m; exp++)
        p *= b;
    
    info->exp = exp;
    info->bb = p;
}
//LINE 1146
static mp_bitcnt_t
mpn_limb_size_in_base_2 (mp_limb_t u)
{
    unsigned shift;
    
    assert (u > 0);
    gmp_clz (shift, u);
    return GMP_LIMB_BITS - shift;
}
//LINE 1156
static size_t
mpn_get_str_bits (unsigned char *sp, unsigned bits, mp_srcptr up, mp_size_t un)
{
    unsigned char mask;
    size_t sn, j;
    mp_size_t i;
    int shift;
    
    sn = ((un - 1) * GMP_LIMB_BITS + mpn_limb_size_in_base_2 (up[un-1])
          + bits - 1) / bits;
    
    mask = (1U << bits) - 1;
    
    for (i = 0, j = sn, shift = 0; j-- > 0;)
    {
        unsigned char digit = up[i] >> shift;
        
        shift += bits;
        
        if (shift >= GMP_LIMB_BITS && ++i < un)
        {
            shift -= GMP_LIMB_BITS;
            digit |= up[i] << (bits - shift);
        }
        sp[j] = digit & mask;
    }
    return sn;
}
//LINE 1187
static size_t
mpn_limb_get_str (unsigned char *sp, mp_limb_t w,
                  const struct gmp_div_inverse *binv)
{
    mp_size_t i;
    for (i = 0; w > 0; i++)
    {
        mp_limb_t h, l, r;
        
        h = w >> (GMP_LIMB_BITS - binv->shift);
        l = w << binv->shift;
        
        gmp_udiv_qrnnd_preinv (w, r, h, l, binv->d1, binv->di);
        assert ( (r << (GMP_LIMB_BITS - binv->shift)) == 0);
        r >>= binv->shift;
        
        sp[i] = r;
    }
    return i;
}
//LINE 1208
static size_t
mpn_get_str_other (unsigned char *sp,
                   int base, const struct mpn_base_info *info,
                   mp_ptr up, mp_size_t un)
{
    struct gmp_div_inverse binv;
    size_t sn;
    size_t i;
    
    mpn_div_qr_1_invert (&binv, base);
    
    sn = 0;
    
    if (un > 1)
    {
        struct gmp_div_inverse bbinv;
        mpn_div_qr_1_invert (&bbinv, info->bb);
        
        do
        {
            mp_limb_t w;
            size_t done;
            w = mpn_div_qr_1_preinv (up, up, un, &bbinv);
            un -= (up[un-1] == 0);
            done = mpn_limb_get_str (sp + sn, w, &binv);
            
            for (sn += done; done < info->exp; done++)
                sp[sn++] = 0;
        }
        while (un > 1);
    }
    sn += mpn_limb_get_str (sp + sn, up[0], &binv);
    
    /* Reverse order */
    for (i = 0; 2*i + 1 < sn; i++)
    {
        unsigned char t = sp[i];
        sp[i] = sp[sn - i - 1];
        sp[sn - i - 1] = t;
    }
    
    return sn;
}
//LINE 1272
static mp_size_t
mpn_set_str_bits (mp_ptr rp, const unsigned char *sp, size_t sn,
                  unsigned bits)
{
    mp_size_t rn;
    size_t j;
    unsigned shift;
    
    for (j = sn, rn = 0, shift = 0; j-- > 0; )
    {
        if (shift == 0)
        {
            rp[rn++] = sp[j];
            shift += bits;
        }
        else
        {
            rp[rn-1] |= (mp_limb_t) sp[j] << shift;
            shift += bits;
            if (shift >= GMP_LIMB_BITS)
            {
                shift -= GMP_LIMB_BITS;
                if (shift > 0)
                    rp[rn++] = (mp_limb_t) sp[j] >> (bits - shift);
            }
        }
    }
    rn = mpn_normalized_size (rp, rn);
    return rn;
}
//LINE 1303
static mp_size_t
mpn_set_str_other (mp_ptr rp, const unsigned char *sp, size_t sn,
                   mp_limb_t b, const struct mpn_base_info *info)
{
    mp_size_t rn;
    mp_limb_t w;
    unsigned k;
    size_t j;
    
    k = 1 + (sn - 1) % info->exp;
    
    j = 0;
    w = sp[j++];
    while (--k != 0)
        w = w * b + sp[j++];
    
    rp[0] = w;
    
    for (rn = (w > 0); j < sn;)
    {
        mp_limb_t cy;
        
        w = sp[j++];
        for (k = 1; k < info->exp; k++)
            w = w * b + sp[j++];
        
        cy = mpn_mul_1 (rp, rp, rn, info->bb);
        cy += mpn_add_1 (rp, rp, rn, w);
        if (cy > 0)
            rp[rn++] = cy;
    }
    assert (j == sn);
    
    return rn;
}
//LINE 1360
/* MPZ interface */
void
mpz_init (mpz_t r)
{
    r->_mp_alloc = 1;
    r->_mp_size = 0;
    r->_mp_d = gmp_xalloc_limbs (1);
}
//LINE 1372
void
mpz_init2 (mpz_t r, mp_bitcnt_t bits)
{
    mp_size_t rn;
    
    bits -= (bits != 0);		/* Round down, except if 0 */
    rn = 1 + bits / GMP_LIMB_BITS;
    
    r->_mp_alloc = rn;
    r->_mp_size = 0;
    r->_mp_d = gmp_xalloc_limbs (rn);
}
//LINE 1384
void
mpz_clear (mpz_t r)
{
    gmp_free (r->_mp_d);
}
//LINE 1390
static mp_ptr
mpz_realloc (mpz_t r, mp_size_t size)
{
    size = GMP_MAX (size, 1);
    
    r->_mp_d = gmp_xrealloc_limbs (r->_mp_d, size);
    r->_mp_alloc = size;
    
    if (GMP_ABS (r->_mp_size) > size)
        r->_mp_size = 0;
    
    return r->_mp_d;
}
//LINE 1404
/* Realloc for an mpz_t WHAT if it has less than NEEDED limbs.  */
#define MPZ_REALLOC(z,n) ((n) > (z)->_mp_alloc			\
? mpz_realloc(z,n)			\
: (z)->_mp_d)
//LINE 1422
void
mpz_set_ui (mpz_t r, unsigned long int x)
{
    if (x > 0)
    {
        r->_mp_size = 1;
        r->_mp_d[0] = x;
    }
    else
        r->_mp_size = 0;
}
//LINE 1410
void
mpz_set_si (mpz_t r, signed long int x)
{
    if (x >= 0)
        mpz_set_ui (r, x);
    else /* (x < 0) */
    {
        r->_mp_size = -1;
        r->_mp_d[0] = GMP_NEG_CAST (unsigned long int, x);
    }
}
//LINE 1434
void
mpz_set (mpz_t r, const mpz_t x)
{
    /* Allow the NOP r == x */
    if (r != x)
    {
        mp_size_t n;
        mp_ptr rp;
        
        n = GMP_ABS (x->_mp_size);
        rp = MPZ_REALLOC (r, n);
        
        mpn_copyi (rp, x->_mp_d, n);
        r->_mp_size = x->_mp_size;
    }
}
//LINE 1465
void
mpz_init_set (mpz_t r, const mpz_t x)
{
    mpz_init (r);
    mpz_set (r, x);
}
//LINE 1713
int
mpz_sgn (const mpz_t u)
{
    mp_size_t usize = u->_mp_size;
    
    return (usize > 0) - (usize < 0);
}
//LINE 1807
void
mpz_swap (mpz_t u, mpz_t v)
{
    MP_SIZE_T_SWAP (u->_mp_size, v->_mp_size);
    MP_SIZE_T_SWAP (u->_mp_alloc, v->_mp_alloc);
    MP_PTR_SWAP (u->_mp_d, v->_mp_d);
}
//LINE 1819
static mp_size_t
mpz_abs_add_ui (mpz_t r, const mpz_t a, unsigned long b)
{
    mp_size_t an;
    mp_ptr rp;
    mp_limb_t cy;
    
    an = GMP_ABS (a->_mp_size);
    if (an == 0)
    {
        r->_mp_d[0] = b;
        return b > 0;
    }
    
    rp = MPZ_REALLOC (r, an + 1);
    
    cy = mpn_add_1 (rp, a->_mp_d, an, b);
    rp[an] = cy;
    an += cy;
    
    return an;
}
//LINE 1844
static mp_size_t
mpz_abs_sub_ui (mpz_t r, const mpz_t a, unsigned long b)
{
    mp_size_t an = GMP_ABS (a->_mp_size);
    mp_ptr rp = MPZ_REALLOC (r, an);
    
    if (an == 0)
    {
        rp[0] = b;
        return -(b > 0);
    }
    else if (an == 1 && a->_mp_d[0] < b)
    {
        rp[0] = b - a->_mp_d[0];
        return -1;
    }
    else
    {
        gmp_assert_nocarry (mpn_sub_1 (rp, a->_mp_d, an, b));
        return mpn_normalized_size (rp, an);
    }
}
//LINE 1867
void
mpz_add_ui (mpz_t r, const mpz_t a, unsigned long b)
{
    if (a->_mp_size >= 0)
        r->_mp_size = mpz_abs_add_ui (r, a, b);
    else
        r->_mp_size = -mpz_abs_sub_ui (r, a, b);
}
//LINE 1876
void
mpz_sub_ui (mpz_t r, const mpz_t a, unsigned long b)
{
    if (a->_mp_size < 0)
        r->_mp_size = -mpz_abs_add_ui (r, a, b);
    else
        r->_mp_size = mpz_abs_sub_ui (r, a, b);
}
//LINE 1894
static mp_size_t
mpz_abs_add (mpz_t r, const mpz_t a, const mpz_t b)
{
    mp_size_t an = GMP_ABS (a->_mp_size);
    mp_size_t bn = GMP_ABS (b->_mp_size);
    mp_ptr rp;
    mp_limb_t cy;
    
    if (an < bn)
    {
        MPZ_SRCPTR_SWAP (a, b);
        MP_SIZE_T_SWAP (an, bn);
    }
    
    rp = MPZ_REALLOC (r, an + 1);
    cy = mpn_add (rp, a->_mp_d, an, b->_mp_d, bn);
    
    rp[an] = cy;
    
    return an + cy;
}
//LINE 1916
static mp_size_t
mpz_abs_sub (mpz_t r, const mpz_t a, const mpz_t b)
{
    mp_size_t an = GMP_ABS (a->_mp_size);
    mp_size_t bn = GMP_ABS (b->_mp_size);
    int cmp;
    mp_ptr rp;
    
    cmp = mpn_cmp4 (a->_mp_d, an, b->_mp_d, bn);
    if (cmp > 0)
    {
        rp = MPZ_REALLOC (r, an);
        gmp_assert_nocarry (mpn_sub (rp, a->_mp_d, an, b->_mp_d, bn));
        return mpn_normalized_size (rp, an);
    }
    else if (cmp < 0)
    {
        rp = MPZ_REALLOC (r, bn);
        gmp_assert_nocarry (mpn_sub (rp, b->_mp_d, bn, a->_mp_d, an));
        return -mpn_normalized_size (rp, bn);
    }
    else
        return 0;
}
//LINE 1941
void
mpz_add (mpz_t r, const mpz_t a, const mpz_t b)
{
    mp_size_t rn;
    
    if ( (a->_mp_size ^ b->_mp_size) >= 0)
        rn = mpz_abs_add (r, a, b);
    else
        rn = mpz_abs_sub (r, a, b);
    
    r->_mp_size = a->_mp_size >= 0 ? rn : - rn;
}
//LINE 1954
void
mpz_sub (mpz_t r, const mpz_t a, const mpz_t b)
{
    mp_size_t rn;
    
    if ( (a->_mp_size ^ b->_mp_size) >= 0)
        rn = mpz_abs_sub (r, a, b);
    else
        rn = mpz_abs_add (r, a, b);
    
    r->_mp_size = a->_mp_size >= 0 ? rn : - rn;
}
//LINE 2006
void
mpz_mul (mpz_t r, const mpz_t u, const mpz_t v)
{
    int sign;
    mp_size_t un, vn, rn;
    mpz_t t;
    mp_ptr tp;
    
    un = u->_mp_size;
    vn = v->_mp_size;
    
    if (un == 0 || vn == 0)
    {
        r->_mp_size = 0;
        return;
    }
    
    sign = (un ^ vn) < 0;
    
    un = GMP_ABS (un);
    vn = GMP_ABS (vn);
    
    mpz_init2 (t, (un + vn) * GMP_LIMB_BITS);
    
    tp = t->_mp_d;
    if (un >= vn)
        mpn_mul (tp, u->_mp_d, un, v->_mp_d, vn);
    else
        mpn_mul (tp, v->_mp_d, vn, u->_mp_d, un);
    
    rn = un + vn;
    rn -= tp[rn-1] == 0;
    
    t->_mp_size = sign ? - rn : rn;
    mpz_swap (r, t);
    mpz_clear (t);
}
//LINE 2120
enum mpz_div_round_mode { GMP_DIV_FLOOR, GMP_DIV_CEIL, GMP_DIV_TRUNC };
//LINE 2123
static int
mpz_div_qr (mpz_t q, mpz_t r,
            const mpz_t n, const mpz_t d, enum mpz_div_round_mode mode)
{
    mp_size_t ns, ds, nn, dn, qs;
    ns = n->_mp_size;
    ds = d->_mp_size;
    
    if (ds == 0)
        gmp_die("mpz_div_qr: Divide by zero.");
    
    if (ns == 0)
    {
        if (q)
            q->_mp_size = 0;
        if (r)
            r->_mp_size = 0;
        return 0;
    }
    
    nn = GMP_ABS (ns);
    dn = GMP_ABS (ds);
    
    qs = ds ^ ns;
    
    if (nn < dn)
    {
        if (mode == GMP_DIV_CEIL && qs >= 0)
        {
            /* q = 1, r = n - d */
            if (r)
                mpz_sub (r, n, d);
            if (q)
                mpz_set_ui (q, 1);
        }
        else if (mode == GMP_DIV_FLOOR && qs < 0)
        {
            /* q = -1, r = n + d */
            if (r)
                mpz_add (r, n, d);
            if (q)
                mpz_set_si (q, -1);
        }
        else
        {
            /* q = 0, r = d */
            if (r)
                mpz_set (r, n);
            if (q)
                q->_mp_size = 0;
        }
        return 1;
    }
    else
    {
        mp_ptr np, qp;
        mp_size_t qn, rn;
        mpz_t tq, tr;
        
        mpz_init_set (tr, n);
        np = tr->_mp_d;
        
        qn = nn - dn + 1;
        
        if (q)
        {
            mpz_init2 (tq, qn * GMP_LIMB_BITS);
            qp = tq->_mp_d;
        }
        else
            qp = NULL;
        
        mpn_div_qr (qp, np, nn, d->_mp_d, dn);
        
        if (qp)
        {
            qn -= (qp[qn-1] == 0);
            
            tq->_mp_size = qs < 0 ? -qn : qn;
        }
        rn = mpn_normalized_size (np, dn);
        tr->_mp_size = ns < 0 ? - rn : rn;
        
        if (mode == GMP_DIV_FLOOR && qs < 0 && rn != 0)
        {
            if (q)
                mpz_sub_ui (tq, tq, 1);
            if (r)
                mpz_add (tr, tr, d);
        }
        else if (mode == GMP_DIV_CEIL && qs >= 0 && rn != 0)
        {
            if (q)
                mpz_add_ui (tq, tq, 1);
            if (r)
                mpz_sub (tr, tr, d);
        }
        
        if (q)
        {
            mpz_swap (tq, q);
            mpz_clear (tq);
        }
        if (r)
            mpz_swap (tr, r);
        
        mpz_clear (tr);
        
        return rn != 0;
    }
}
//LINE 2247
void
mpz_tdiv_qr (mpz_t q, mpz_t r, const mpz_t n, const mpz_t d)
{
    mpz_div_qr (q, r, n, d, GMP_DIV_TRUNC);
}
//LINE 3977
/* MPZ base conversion. */

size_t
mpz_sizeinbase (const mpz_t u, int base)
{
    mp_size_t un;
    mp_srcptr up;
    mp_ptr tp;
    mp_bitcnt_t bits;
    struct gmp_div_inverse bi;
    size_t ndigits;
    
    assert (base >= 2);
    assert (base <= 36);
    
    un = GMP_ABS (u->_mp_size);
    if (un == 0)
        return 1;
    
    up = u->_mp_d;
    
    bits = (un - 1) * GMP_LIMB_BITS + mpn_limb_size_in_base_2 (up[un-1]);
    switch (base)
    {
        case 2:
            return bits;
        case 4:
            return (bits + 1) / 2;
        case 8:
            return (bits + 2) / 3;
        case 16:
            return (bits + 3) / 4;
        case 32:
            return (bits + 4) / 5;
            /* FIXME: Do something more clever for the common case of base
             10. */
    }
    
    tp = gmp_xalloc_limbs (un);
    mpn_copyi (tp, up, un);
    mpn_div_qr_1_invert (&bi, base);
    
    ndigits = 0;
    do
    {
        ndigits++;
        mpn_div_qr_1_preinv (tp, tp, un, &bi);
        un -= (tp[un-1] == 0);
    }
    while (un > 0);
    
    gmp_free (tp);
    return ndigits;
}
//LINE 4032
char *
mpz_get_str (char *sp, int base, const mpz_t u)
{
    unsigned bits;
    const char *digits;
    mp_size_t un;
    size_t i, sn;
    
    if (base >= 0)
    {
        digits = "0123456789abcdefghijklmnopqrstuvwxyz";
    }
    else
    {
        base = -base;
        digits = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
    }
    if (base <= 1)
        base = 10;
    if (base > 36)
        return NULL;
    
    sn = 1 + mpz_sizeinbase (u, base);
    if (!sp)
        sp = (char *) gmp_xalloc (1 + sn);
    
    un = GMP_ABS (u->_mp_size);
    
    if (un == 0)
    {
        sp[0] = '0';
        sp[1] = '\0';
        return sp;
    }
    
    i = 0;
    
    if (u->_mp_size < 0)
        sp[i++] = '-';
    
    bits = mpn_base_power_of_two_p (base);
    
    if (bits)
    /* Not modified in this case. */
        sn = i + mpn_get_str_bits ((unsigned char *) sp + i, bits, u->_mp_d, un);
    else
    {
        struct mpn_base_info info;
        mp_ptr tp;
        
        mpn_get_base_info (&info, base);
        tp = gmp_xalloc_limbs (un);
        mpn_copyi (tp, u->_mp_d, un);
        
        sn = i + mpn_get_str_other ((unsigned char *) sp + i, base, &info, tp, un);
        gmp_free (tp);
    }
    
    for (; i < sn; i++)
        sp[i] = digits[(unsigned char) sp[i]];
    
    sp[sn] = '\0';
    return sp;
}
//LINE 4097
int
mpz_set_str (mpz_t r, const char *sp, int base)
{
    unsigned bits;
    mp_size_t rn, alloc;
    mp_ptr rp;
    size_t sn;
    int sign;
    unsigned char *dp;
    
    assert (base == 0 || (base >= 2 && base <= 36));
    
    while (isspace( (unsigned char) *sp))
        sp++;
    
    sign = (*sp == '-');
    sp += sign;
    
    if (base == 0)
    {
        if (*sp == '0')
        {
            sp++;
            if (*sp == 'x' || *sp == 'X')
            {
                base = 16;
                sp++;
            }
            else if (*sp == 'b' || *sp == 'B')
            {
                base = 2;
                sp++;
            }
            else
                base = 8;
        }
        else
            base = 10;
    }
    
    sn = strlen (sp);
    dp = (unsigned char *) gmp_xalloc (sn + (sn == 0));
    
    for (sn = 0; *sp; sp++)
    {
        unsigned digit;
        
        if (isspace ((unsigned char) *sp))
            continue;
        if (*sp >= '0' && *sp <= '9')
            digit = *sp - '0';
        else if (*sp >= 'a' && *sp <= 'z')
            digit = *sp - 'a' + 10;
        else if (*sp >= 'A' && *sp <= 'Z')
            digit = *sp - 'A' + 10;
        else
            digit = base; /* fail */
        
        if (digit >= base)
        {
            gmp_free (dp);
            r->_mp_size = 0;
            return -1;
        }
        
        dp[sn++] = digit;
    }
    
    bits = mpn_base_power_of_two_p (base);
    
    if (bits > 0)
    {
        alloc = (sn * bits + GMP_LIMB_BITS - 1) / GMP_LIMB_BITS;
        rp = MPZ_REALLOC (r, alloc);
        rn = mpn_set_str_bits (rp, dp, sn, bits);
    }
    else
    {
        struct mpn_base_info info;
        mpn_get_base_info (&info, base);
        alloc = (sn + info.exp - 1) / info.exp;
        rp = MPZ_REALLOC (r, alloc);
        rn = mpn_set_str_other (rp, dp, sn, base, &info);
    }
    assert (rn <= alloc);
    gmp_free (dp);
    
    r->_mp_size = sign ? - rn : rn;
    
    return 0;
}
//LINE 4189
int
mpz_init_set_str (mpz_t r, const char *sp, int base)
{
    mpz_init (r);
    return mpz_set_str (r, sp, base);
}
/////////////////
#include <stdlib.h>

#define OP_ADD 1
#define OP_SUB 2
#define OP_MUL 4
#define OP_DIV 8
FILE *fp;
FILE *fpo;
char temp1[1000000] = {0};
char temp2[1000000] = {0};
char temp3[1000000] = {0};
char lhs_str[1000000] = {0};
char rhs_str[1000000] = {0};
char ans_str[1000000] = {0};
char ans1_str[1000000] = {0};
void error(void){
    fprintf(fpo,"Error input\n");
}
int parseOp(void){
    int retval;
    char op;
    sscanf(temp3,"%c",&op);
    switch (op) {
        case '+':
            retval = OP_ADD;
            break;
        case '-':
            retval = OP_SUB;
            break;
        case '*':
            retval = OP_MUL;
            break;
        case '/':
            retval = OP_DIV;
            break;
        default:
            retval = -1;
            break;
    }
    return retval;
}
int parseCmdArg(int argc,const char * argv[]){
    if (argc != 3) {
        printf("please input a file name as the input file such as \"input.txt\"\n");
        return -1;
    }
    return 0;
}
int u_fscanfs(FILE *,char temp[]){
    long long i = 0;
    char u_temp;
    //step1 get first character
    u_temp = fgetc(fp);
    if (u_temp == EOF) {
        return EOF;
    }
    if (u_temp != '\n') {
        temp[i] = u_temp;
        i++;
    }else{
        temp[i] = '\0';
        return 0;
    }
    //step2 get others
    while ((u_temp = fgetc(fp)) != '\n' && u_temp != EOF) {
        temp[i] = u_temp;
        i++;
    }
    temp[i] = '\0';
    return 1;
}
int main(int argc, const char * argv[]) {
    mpz_t lhs,rhs,ans,reminder;
    int retval = 1;
    int retL = 1;
    int retR = 1;
    mpz_init(lhs);
    mpz_init(rhs);
    mpz_init(ans);
    mpz_init(reminder);
    retval = parseCmdArg(argc, argv);
    if (retval != 0) {
        return -1;
    }
    fp = fopen(argv[1],"r");
    fpo = fopen(argv[2], "w");
    if (!fp) {
        printf("con't open the input file\n");
        return -1;
    }
    while (1) {
        retval = u_fscanfs(fp,temp1);
        if (retval == EOF) {
            break;
        }
        retL = u_fscanfs(fp,temp2);
        retR = u_fscanfs(fp,temp3);
        if(retL != 1 || retR != 1){
            //just to handle empty line
            error();
            continue;
        }
        // to handel space in the middle of line
        char has_space_worng_L = 0;
        char has_space_worng_R = 0;
        retL = sscanf(temp1,"%s%c",lhs_str,&has_space_worng_L);
        retR = sscanf(temp2,"%s%c",rhs_str,&has_space_worng_R);
        
        if(retL != 1 || retR != 1){
            error();
            continue;
        }
        retL = mpz_init_set_str(lhs, lhs_str, 10);
        retR = mpz_init_set_str(rhs, rhs_str, 10);
        if (retL == -1 || retR == -1) {
            error();
            continue;
        }
        
        retval = parseOp();
        if (retval == -1) {
            error();
            continue;
        }
        
        if (retval == OP_SUB) {
            mpz_sub(ans,lhs,rhs);
            fprintf(fpo,"%s\n",mpz_get_str(ans_str,10, ans));
        }
        else if (retval == OP_ADD) {
            mpz_add(ans,lhs,rhs);
            fprintf(fpo,"%s\n",mpz_get_str(ans_str,10, ans));
        }
        else if (retval == OP_MUL) {
            mpz_mul(ans,lhs,rhs);
            fprintf(fpo,"%s\n",mpz_get_str(ans_str,10, ans));
        }
        else if (retval == OP_DIV) {
            if (!mpz_sgn(rhs)) {
                //if divisor is zero
                error();
                continue;
            }
            if(mpz_sgn(rhs) == -1 || mpz_sgn(lhs) == -1){
                error();
                continue;
            }
            mpz_tdiv_qr(ans,reminder,lhs,rhs);
            fprintf(fpo,"%s %s\n",mpz_get_str(ans_str,10, ans),mpz_get_str(ans1_str, 10,reminder));
        }
    }
    mpz_clear(ans);
    mpz_clear(lhs);
    mpz_clear(rhs);
    mpz_clear(reminder);
    fclose(fp);
    fclose(fpo);
    return  0;
}
