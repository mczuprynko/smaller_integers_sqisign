#include "hnf_internal.h"
#include "internal.h"
#include <stdio.h>
#include "bench.h"

// HNF test function
int
ibz_mat_4x4_is_hnf(const ibz_mat_4x4_t *mat)
{
    int res = 1;
    int found;
    int ind = 0;
    ibz_t zero;
    ibz_init(&zero);
    // upper triangular
    for (int i = 0; i < 4; i++) {
        // upper triangular
        for (int j = 0; j < i; j++) {
            res = res && ibz_is_zero(&((*mat)[i][j]));
        }
        // find first non 0 element of line
        found = 0;
        for (int j = i; j < 4; j++) {
            if (found) {
                // all values are positive, and first non-0 is the largest of that line
                res = res && (ibz_cmp(&((*mat)[i][j]), &zero) >= 0);
                res = res && (ibz_cmp(&((*mat)[i][ind]), &((*mat)[i][j])) > 0);
            } else {
                if (!ibz_is_zero(&((*mat)[i][j]))) {
                    found = 1;
                    ind = j;
                    // mustbe non-negative
                    res = res && (ibz_cmp(&((*mat)[i][j]), &zero) > 0);
                }
            }
        }
    }
    // check that first nom-zero elements ndex per column is strictly increasing
    int linestart = -1;
    int i = 0;
    for (int j = 0; j < 4; j++) {
        while ((i < 4) && (ibz_is_zero(&((*mat)[i][j])))) {
            i = i + 1;
        }
        if (i != 4) {
            res = res && (linestart < i);
        }
        i = 0;
    }
    ibz_finalize(&zero);
    return res;
}

// Untested HNF helpers
//  centered mod
void
ibz_vec_4_linear_combination_mod(ibz_vec_4_t *lc,
                                 const ibz_t *coeff_a,
                                 const ibz_vec_4_t *vec_a,
                                 const ibz_t *coeff_b,
                                 const ibz_vec_4_t *vec_b,
                                 const ibz_t *mod)
{
    ibz_t prod, m;
    ibz_vec_4_t sums;
    ibz_vec_4_init(&sums);
    ibz_init(&prod);
    ibz_init(&m);
    ibz_copy(&m, mod);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&(sums[i]), coeff_a, &((*vec_a)[i]));
        ibz_mul(&prod, coeff_b, &((*vec_b)[i]));
        ibz_add(&(sums[i]), &(sums[i]), &prod);
        ibz_centered_mod(&(sums[i]), &(sums[i]), &m);
    }
    for (int i = 0; i < 4; i++) {
        ibz_copy(&((*lc)[i]), &(sums[i]));
    }
    ibz_finalize(&prod);
    ibz_finalize(&m);
    ibz_vec_4_finalize(&sums);
}

void
ibz_vec_4_copy_mod(ibz_vec_4_t *res, const ibz_vec_4_t *vec, const ibz_t *mod)
{
    ibz_t m;
    ibz_init(&m);
    ibz_copy(&m, mod);
    for (int i = 0; i < 4; i++) {
        ibz_centered_mod(&((*res)[i]), &((*vec)[i]), &m);
    }
    ibz_finalize(&m);
}

// no need to center this, and not 0
void
ibz_vec_4_scalar_mul_mod(ibz_vec_4_t *prod, const ibz_t *scalar, const ibz_vec_4_t *vec, const ibz_t *mod)
{
    ibz_t m, s;
    ibz_init(&m);
    ibz_init(&s);
    ibz_copy(&s, scalar);
    ibz_copy(&m, mod);
    for (int i = 0; i < 4; i++) {
        ibz_mul(&((*prod)[i]), &((*vec)[i]), &s);
        ibz_mod(&((*prod)[i]), &((*prod)[i]), &m);
    }
    ibz_finalize(&m);
    ibz_finalize(&s);
}

// Algorithm used is the one at number 2.4.8 in Henri Cohen's "A Course in Computational Algebraic
// Number Theory" (Springer Verlag, in series "Graduate texts in Mathematics") from 1993
//  assumes ibz_xgcd outputs u,v which are small in absolute value (as described in the
//  book)
void
ibz_mat_4xn_hnf_mod_core(ibz_mat_4x4_t *hnf, int generator_number, const ibz_vec_4_t *generators, const ibz_t *mod)
{
    int i = 3;
    assert(generator_number > 3);
    int n = generator_number;
    int j = n - 1;
    int k = n - 1;
    ibz_t b, u, v, d, q, m, coeff_1, coeff_2, r;
    ibz_vec_4_t c;
    ibz_vec_4_t a[generator_number];
    ibz_vec_4_t w[4];
    ibz_init(&b);
    ibz_init(&d);
    ibz_init(&u);
    ibz_init(&v);
    ibz_init(&r);
    ibz_init(&m);
    ibz_init(&q);
    ibz_init(&coeff_1);
    ibz_init(&coeff_2);
    ibz_vec_4_init(&c);
    PRINT_BACKTRACE_SYMBOLS(0LU);
    for (int h = 0; h < n; h++) {
        if (h < 4)
            ibz_vec_4_init(&(w[h]));
        ibz_vec_4_init(&(a[h]));
        ibz_copy(&(a[h][0]), &(generators[h][0]));
        ibz_copy(&(a[h][1]), &(generators[h][1]));
        ibz_copy(&(a[h][2]), &(generators[h][2]));
        ibz_copy(&(a[h][3]), &(generators[h][3]));
    }
    assert(ibz_cmp(mod, &ibz_const_zero) > 0);
    ibz_copy(&m, mod);
    while (i != -1) {
        while (j != 0) {
            j = j - 1;
            if (!ibz_is_zero(&(a[j][i]))) {
                // assumtion that ibz_xgcd outputs u,v which are small in absolute
                // value is needed here also, needs u non 0, but v can be 0 if needed
                ibz_xgcd_with_u_not_0(&d, &u, &v, &(a[k][i]), &(a[j][i]));
                ibz_vec_4_linear_combination(&c, &u, &(a[k]), &v, &(a[j]));
                ibz_div(&coeff_1, &r, &(a[k][i]), &d);
                ibz_div(&coeff_2, &r, &(a[j][i]), &d);
                ibz_neg(&coeff_2, &coeff_2);
                ibz_vec_4_linear_combination_mod(
                    &(a[j]), &coeff_1, &(a[j]), &coeff_2, &(a[k]), &m); // do lin comb mod m
                ibz_vec_4_copy_mod(&(a[k]), &c, &m);                    // mod m in copy
            }
        }
        ibz_xgcd_with_u_not_0(&d, &u, &v, &(a[k][i]), &m);
        ibz_vec_4_scalar_mul_mod(&(w[i]), &u, &(a[k]), &m); // mod m in scalar mult
        if (ibz_is_zero(&(w[i][i]))) {
            ibz_copy(&(w[i][i]), &m);
        }
        for (int h = i + 1; h < 4; h++) {
            ibz_div_floor(&q, &r, &(w[h][i]), &(w[i][i]));
            ibz_neg(&q, &q);
            ibz_vec_4_linear_combination(&(w[h]), &ibz_const_one, &(w[h]), &q, &(w[i]));
        }
        ibz_div(&m, &r, &m, &d);
        assert(ibz_is_zero(&r));
        if (i != 0) {
            k = k - 1;
            i = i - 1;
            j = k;
            if (ibz_is_zero(&(a[k][i])))
                ibz_copy(&(a[k][i]), &m);

        } else {
            k = k - 1;
            i = i - 1;
            j = k;
        }
    }
    for (j = 0; j < 4; j++) {
        for (i = 0; i < 4; i++) {
            ibz_copy(&((*hnf)[i][j]), &(w[j][i]));
        }
    }

    ibz_finalize(&b);
    ibz_finalize(&d);
    ibz_finalize(&u);
    ibz_finalize(&v);
    ibz_finalize(&r);
    ibz_finalize(&q);
    ibz_finalize(&coeff_1);
    ibz_finalize(&coeff_2);
    ibz_finalize(&m);
    ibz_vec_4_finalize(&c);
    for (int h = 0; h < n; h++) {
        if (h < 4)
            ibz_vec_4_finalize(&(w[h]));
        ibz_vec_4_finalize(&(a[h]));
    }
}

void
ibz_mat_4xn_hnf_simple_mod_core(ibz_mat_4x4_t *hnf, int generator_number, const ibz_vec_4_t *generators, const ibz_t *mod)
{
    assert(generator_number > 3);
    assert(!ibz_is_zero(mod));
    ibz_t tmp, d, u, v, s, t;
    ibz_vec_4_t a[generator_number + 4], tmp_vec;
    ibz_init(&tmp); ibz_init(&d); ibz_init(&u); ibz_init(&v); ibz_init(&s); ibz_init(&t);
    ibz_vec_4_init(&tmp_vec);
    for (int i = 0; i < generator_number; i++) {
        ibz_vec_4_init(&a[i]);
        ibz_vec_4_copy_mod(&a[i], &generators[i], mod);
    }
    for (int i = generator_number; i < generator_number + 4; i++){
        ibz_vec_4_init(&a[i]);
        ibz_copy(&a[i][3 - (i - generator_number)], mod);
    }

    for (int i = 3; i >= 0; i--){
        for (int j = 4-i; j < generator_number + 4 - i; j++) {
            if (ibz_is_zero(&(a[j][i]))) continue;
            
            ibz_xgcd_with_u_not_0(&d, &u, &v, &a[3-i][i], &a[j][i]);
            ibz_div(&s, &tmp, &(a[3-i][i]), &d);
            assert(ibz_is_zero(&tmp));
            ibz_div(&t, &tmp, &(a[j][i]), &d);
            assert(ibz_is_zero(&tmp));
            ibz_neg(&s, &s);
            ibz_vec_4_linear_combination(&tmp_vec, &u, &a[3-i], &v, &a[j]);
            ibz_vec_4_linear_combination_mod(&a[j], &t, &a[3-i], &s, &a[j], mod);
            ibz_copy(&a[3-i][i], &tmp_vec[i]);
            for (int k = 0; k < i; k++){
                ibz_mod(&a[3-i][k], &tmp_vec[k], mod);
            }
        }
        

        for (int j = 0; j < 3-i; j++) {
            ibz_div_floor(&d, &a[j][i], &a[j][i], &a[3-i][i]);
            for (int k = 0; k < i; k++){
                ibz_mul(&tmp, &d, &a[3-i][k]);
                ibz_sub(&a[j][k], &a[j][k], &tmp);
                ibz_mod(&a[j][k], &a[j][k], mod);
            }
        }
    }

    for (int i = 4; i < generator_number + 4; i++){
        assert(ibz_vec_4_is_zero(&a[i]));
    }

    for (int i = 0; i < 4; i++){
        for (int j = 0; j < 4; j++){
            ibz_copy(&(*hnf)[j][3-i], &a[i][j]);
        }
    }

    ibz_finalize(&tmp); ibz_finalize(&d); ibz_finalize(&u); ibz_finalize(&v); ibz_finalize(&s); ibz_finalize(&t);
    ibz_vec_4_finalize(&tmp_vec);
    for (int i = 0; i < generator_number + 4; i++){
        ibz_vec_4_finalize(&a[i]);
    }
}