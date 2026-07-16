#include "internal.h"
#include "intbig.h"
#include <stdio.h>

void cbz_add(cbz_t *c, const cbz_t *a, const cbz_t *b){
    ibz_add(&(*c)[0], &(*a)[0], &(*b)[0]);
    ibz_add(&(*c)[1], &(*a)[1], &(*b)[1]);
}

void cbz_sub(cbz_t *c, const cbz_t *a, const cbz_t *b){
    ibz_sub(&(*c)[0], &(*a)[0], &(*b)[0]);
    ibz_sub(&(*c)[1], &(*a)[1], &(*b)[1]);
}

void cbz_mult(cbz_t *c, const cbz_t *a, const cbz_t *b){
    ibz_t prod00, prod01, prod10, prod11;
    
    ibz_init(&prod00); ibz_init(&prod01); ibz_init(&prod10); ibz_init(&prod11);
    ibz_mul(&prod00, &(*a)[0], &(*b)[0]);
    ibz_mul(&prod01, &(*a)[0], &(*b)[1]);
    ibz_mul(&prod10, &(*a)[1], &(*b)[0]);
    ibz_mul(&prod11, &(*a)[1], &(*b)[1]);
    ibz_sub(&(*c)[0], &prod00, &prod11);
    ibz_add(&(*c)[1], &prod01, &prod10);

    ibz_finalize(&prod00); ibz_finalize(&prod01); ibz_finalize(&prod10); ibz_finalize(&prod11);
}

void cbz_swap(cbz_t *a, cbz_t *b){
    ibz_swap(&(*a)[0], &(*b)[0]);
    ibz_swap(&(*a)[1], &(*b)[1]);
}

void cbz_vec_2_swap(cbz_vec_2_t *a, cbz_vec_2_t *b){
    cbz_swap(&(*a)[0], &(*b)[0]);
    cbz_swap(&(*a)[1], &(*b)[1]);
}

void cbz_set_uint32(cbz_t *a, const int32_t rel, const int32_t im) {
    ibz_set(&(*a)[0], rel);
    ibz_set(&(*a)[1], im);
}

void cbz_set_ibz(cbz_t *a, const ibz_t *rel, const ibz_t *im) {
    ibz_copy(&(*a)[0], rel);
    ibz_copy(&(*a)[1], im);
}

void cbz_conjugate(cbz_t *r, const cbz_t *a){
    ibz_copy(&(*r)[0], &(*a)[0]);
    ibz_neg(&(*r)[1], &(*a)[1]);
}

void
ibz_rounded_div(ibz_t *q, ibz_t *r, const ibz_t *a, const ibz_t *b)
{
    ibz_t half_b, tmp;
    ibz_init(&half_b); ibz_init(&tmp);

    mpz_div_2exp(half_b, *b, 1);
    ibz_add(&tmp, a, &half_b);
    ibz_div_floor(q, r, &tmp, b);
    ibz_sub(r, r, &half_b);

    ibz_finalize(&half_b); ibz_finalize(&tmp);
}

void cbz_rounded_div_ibz(cbz_t *q, cbz_t *r, const cbz_t *a, const ibz_t *div){
    ibz_rounded_div(&(*q)[0], &(*r)[0], &(*a)[0], div);
    ibz_rounded_div(&(*q)[1], &(*r)[1], &(*a)[1], div);
}



void cbz_init(cbz_t *x){
    ibz_init(&(*x)[0]);
    ibz_init(&(*x)[1]);
}

void cbz_vec_2_init(cbz_vec_2_t *x){
    cbz_init(&(*x)[0]);
    cbz_init(&(*x)[1]);
}

void cbz_mat_2x2_init(cbz_mat_2x2_t *x){
    cbz_vec_2_init(&(*x)[0]);
    cbz_vec_2_init(&(*x)[1]);
}

void cbz_finalize(cbz_t *x){
    ibz_finalize(&(*x)[0]);
    ibz_finalize(&(*x)[1]);
}

void cbz_vec_2_finalize(cbz_vec_2_t *x){
    cbz_finalize(&(*x)[0]);
    cbz_finalize(&(*x)[1]);
}

void cbz_mat_2x2_finalize(cbz_mat_2x2_t *x){
    cbz_vec_2_finalize(&(*x)[0]);
    cbz_vec_2_finalize(&(*x)[1]);
}


void inter_gauss_lagrange_4x2_swap(cbz_mat_2x2_t *uni_mat, cbz_mat_2x2_t *uv, cbz_mat_2x2_t *gram){
    cbz_swap(&(*gram)[0][0], &(*gram)[1][1]);
    cbz_swap(&(*gram)[1][0], &(*gram)[0][1]);

    cbz_swap(&(*uni_mat)[0][0], &(*uni_mat)[1][0]);
    cbz_swap(&(*uni_mat)[0][1], &(*uni_mat)[1][1]);

    cbz_swap(&(*uv)[0][0], &(*uv)[1][0]);
    cbz_swap(&(*uv)[0][1], &(*uv)[1][1]);
}

void cbz_copy(cbz_t *r, const cbz_t *a){
    ibz_copy(&(*r)[0], &(*a)[0]);
    ibz_copy(&(*r)[1], &(*a)[1]);
}

void cbz_neg(cbz_t *r, const cbz_t *a){
    ibz_neg(&(*r)[0], &(*a)[0]);
    ibz_neg(&(*r)[1], &(*a)[1]);
}

void cbz_vec_2_copy(cbz_vec_2_t *r, const cbz_vec_2_t *a){
    cbz_copy(&(*r)[0], &(*a)[0]);
    cbz_copy(&(*r)[1], &(*a)[1]);
}

void cbz_copy_quat(cbz_vec_2_t *res, const ibz_vec_4_t *x) {
    ibz_copy(&(*res)[0][0], &(*x)[0]);
    ibz_copy(&(*res)[0][1], &(*x)[1]);
    ibz_copy(&(*res)[1][0], &(*x)[2]);
    ibz_copy(&(*res)[1][1], &(*x)[3]);
}

void cbz_rounded_div(cbz_t *q, cbz_t *r, const cbz_t *a, const cbz_t *div){
    cbz_t c_div, norm, tmp, tmp_r;
    cbz_init(&c_div); cbz_init(&norm); cbz_init(&tmp); cbz_init(&tmp_r);

    cbz_conjugate(&c_div, div);
    cbz_mult(&norm, &c_div, div);
    cbz_mult(&c_div, a, &c_div);
    cbz_rounded_div_ibz(&tmp, &tmp_r, &c_div, &norm[0]);
    cbz_mult(&norm, div, &tmp);
    cbz_sub(r, a, &norm);
    cbz_copy(q, &tmp);

    cbz_finalize(&c_div); cbz_finalize(&norm); cbz_finalize(&tmp); cbz_finalize(&tmp_r);
}

void cbz_set_quat(ibz_vec_4_t *res, const cbz_vec_2_t *x) {
    ibz_copy(&(*res)[0], &(*x)[0][0]);
    ibz_copy(&(*res)[1], &(*x)[0][1]);
    ibz_copy(&(*res)[2], &(*x)[1][0]);
    ibz_copy(&(*res)[3], &(*x)[1][1]);
}

void cbz_vec_2_conjugate(cbz_vec_2_t *r, const cbz_vec_2_t *a){
    cbz_conjugate(&(*r)[0], &(*a)[0]);
    cbz_conjugate(&(*r)[1], &(*a)[1]);
}

void cbz_vec_2_mult_i(cbz_vec_2_t *r, const cbz_vec_2_t *a){
    ibz_neg(&(*r)[0][0], &(*a)[0][1]);
    ibz_copy(&(*r)[0][1], &(*a)[0][0]);
    ibz_neg(&(*r)[1][0], &(*a)[1][1]);
    ibz_copy(&(*r)[1][1], &(*a)[1][0]);
}

void cbz_mat_2x2_from_hnf_ibz_lat(cbz_mat_2x2_t *uv, const ibz_mat_4x4_t *mat){

    cbz_vec_2_t *u = &(*uv)[0];
    cbz_vec_2_t *v = &(*uv)[1];
    
    ibz_copy(&(*u)[0][0], &(*mat)[0][3]);
    ibz_copy(&(*u)[0][1], &(*mat)[1][3]);
    ibz_copy(&(*u)[1][0], &(*mat)[2][3]);
    ibz_copy(&(*u)[1][1], &(*mat)[3][3]);


    ibz_copy(&(*v)[0][0], &(*mat)[0][1]);
    ibz_copy(&(*v)[0][1], &(*mat)[1][1]);
    ibz_copy(&(*v)[1][0], &(*mat)[2][1]);
    ibz_copy(&(*v)[1][1], &(*mat)[3][1]);
}

void ibz_lat_from_cbz_mat_2x2(ibz_mat_4x4_t *mat, const cbz_mat_2x2_t *uv){

    cbz_vec_2_t i_u, i_v;
    const cbz_vec_2_t *u = &(*uv)[0];
    const cbz_vec_2_t *v = &(*uv)[1];

    cbz_vec_2_init(&i_u); cbz_vec_2_init(&i_v);

    cbz_vec_2_mult_i(&i_u, u);
    cbz_vec_2_mult_i(&i_v, v);
    
    ibz_copy(&(*mat)[0][0], &(*u)[0][0]);
    ibz_copy(&(*mat)[1][0], &(*u)[0][1]);
    ibz_copy(&(*mat)[2][0], &(*u)[1][0]);
    ibz_copy(&(*mat)[3][0], &(*u)[1][1]);

    ibz_copy(&(*mat)[0][1], &i_u[0][0]);
    ibz_copy(&(*mat)[1][1], &i_u[0][1]);
    ibz_copy(&(*mat)[2][1], &i_u[1][0]);
    ibz_copy(&(*mat)[3][1], &i_u[1][1]);

    ibz_copy(&(*mat)[0][2], &(*v)[0][0]);
    ibz_copy(&(*mat)[1][2], &(*v)[0][1]);
    ibz_copy(&(*mat)[2][2], &(*v)[1][0]);
    ibz_copy(&(*mat)[3][2], &(*v)[1][1]);

    ibz_copy(&(*mat)[0][3], &i_v[0][0]);
    ibz_copy(&(*mat)[1][3], &i_v[0][1]);
    ibz_copy(&(*mat)[2][3], &i_v[1][0]);
    ibz_copy(&(*mat)[3][3], &i_v[1][1]);

    cbz_vec_2_finalize(&i_u); cbz_vec_2_finalize(&i_v);
}

void inter_gauss_lagrange_4x2_reduce(cbz_mat_2x2_t *uni_mat, cbz_mat_2x2_t *uv, cbz_mat_2x2_t *gram){
    cbz_t conj_fact, fact, prod0, prod1, tmp;
    cbz_vec_2_t *u = &(*uv)[0];
    cbz_vec_2_t *v = &(*uv)[1];

    cbz_init(&fact); cbz_init(&prod0); cbz_init(&prod1); cbz_init(&conj_fact); cbz_init(&tmp);
    
    // compute f = round(<u, v>/ |u|^2)
    if (!ibz_is_zero(&(*gram)[0][0][0])) {
        cbz_rounded_div_ibz(&fact, &tmp, &(*gram)[0][1], &(*gram)[0][0][0]);
    }
    cbz_conjugate(&fact, &fact);

    // compute v = v - f*u
    cbz_mult(&prod0, &fact, &(*u)[0]);
    cbz_mult(&prod1, &fact, &(*u)[1]);
    cbz_sub(&(*v)[0], &(*v)[0], &prod0);
    cbz_sub(&(*v)[1], &(*v)[1], &prod1);

    // update |v|^2 = |v|^2 + |f|^2 * |u|^2 - 2*Re(f*<u, v>)
    cbz_conjugate(&conj_fact, &fact);
    cbz_mult(&prod0, &conj_fact, &fact);
    cbz_mult(&prod0, &prod0, &(*gram)[0][0]);
    cbz_add(&(*gram)[1][1], &(*gram)[1][1], &prod0);

    cbz_mult(&prod0, &fact, &(*gram)[0][1]);
    mpz_mul_2exp(prod0[0], prod0[0], 1);
    ibz_sub(&(*gram)[1][1][0], &(*gram)[1][1][0], &prod0[0]);

    // update <u, v>  = <u, v> - conj(f) * |u|
    cbz_mult(&prod0, &(*gram)[0][0], &conj_fact);
    cbz_sub(&(*gram)[0][1], &(*gram)[0][1], &prod0);

    // update <v, u>  = conj(<u, v>)
    cbz_conjugate(&(*gram)[1][0], &(*gram)[0][1]);

    // update uni_mat by [[1, 0], [-f, 1]]
    cbz_mult(&prod0, &fact, &(*uni_mat)[0][0]);
    cbz_mult(&prod1, &fact, &(*uni_mat)[0][1]);
    cbz_sub(&(*uni_mat)[1][0], &(*uni_mat)[1][0], &prod0);
    cbz_sub(&(*uni_mat)[1][1], &(*uni_mat)[1][1], &prod1);

    cbz_finalize(&fact); cbz_finalize(&prod0); cbz_finalize(&prod1); cbz_finalize(&conj_fact); cbz_finalize(&tmp);
}

void cbz_mult_ibz(cbz_t *c, const cbz_t *a, const ibz_t *b){
    ibz_mul(&(*c)[0], &(*a)[0], b);
    ibz_mul(&(*c)[1], &(*a)[1], b);
}

void cbz_div_ibz(cbz_t *c, cbz_t *r, const cbz_t *a, const ibz_t *b){
    ibz_div(&(*c)[0], &(*r)[0], &(*a)[0], b);
    ibz_div(&(*c)[1], &(*r)[1], &(*a)[1], b);
}

void 
cbz_div_floor_norm(cbz_t *c, cbz_t *r, const cbz_t *a, const cbz_t *b) {
    cbz_t conj_b, norm_b, prod;
    cbz_init(&conj_b); cbz_init(&norm_b); cbz_init(&prod);

    cbz_conjugate(&conj_b, b);
    cbz_mult(&norm_b, &conj_b, b);
    cbz_mult(&prod, a, &conj_b);
    cbz_div_ibz(c, r, &prod, &norm_b[0]);

    cbz_finalize(&conj_b); cbz_finalize(&norm_b);  cbz_finalize(&prod);
}

void
cbz_lideal_class_gram(cbz_mat_2x2_t *uv, cbz_mat_2x2_t *G, const quat_left_ideal_t *lideal, const quat_alg_t *alg)
{
    cbz_mat_2x2_from_hnf_ibz_lat(uv, &lideal->lattice.basis);
    cbz_gram_2x2(G, uv, alg);

    // divide by norm · denominator²
    ibz_t divisor;
    cbz_t rmd;
    ibz_init(&divisor);
    cbz_init(&rmd);

    ibz_mul(&divisor, &(lideal->lattice.denom), &(lideal->lattice.denom));
    ibz_mul(&divisor, &divisor, &(lideal->norm));
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            cbz_div_ibz(&(*G)[i][j], &rmd, &(*G)[i][j], &divisor);
            assert(ibz_is_zero(&rmd[0]));
            assert(ibz_is_zero(&rmd[1]));
        }
    }

    cbz_finalize(&rmd);
    ibz_finalize(&divisor);
}

void cbz_gram_2x2(cbz_mat_2x2_t *gram, cbz_mat_2x2_t *uv, const quat_alg_t *alg){
    cbz_t conj0, conj1, prod;
    cbz_vec_2_t *u = &(*uv)[0];
    cbz_vec_2_t *v = &(*uv)[1];

    cbz_init(&conj0); cbz_init(&conj1); cbz_init(&prod);

    cbz_conjugate(&conj0, &(*u)[0]);
    cbz_conjugate(&conj1, &(*u)[1]);
    cbz_mult(&(*gram)[0][0], &(*u)[0], &conj0);
    cbz_mult(&prod, &(*u)[1], &conj1);
    cbz_mult_ibz(&prod, &prod, &(alg->p));
    cbz_add(&(*gram)[0][0], &(*gram)[0][0], &prod);

    cbz_mult(&(*gram)[1][0], &(*v)[0], &conj0);
    cbz_mult(&prod, &(*v)[1], &conj1);
    cbz_mult_ibz(&prod, &prod, &(alg->p));
    cbz_add(&(*gram)[1][0], &(*gram)[1][0], &prod);

    cbz_conjugate(&conj0, &(*v)[0]);
    cbz_conjugate(&conj1, &(*v)[1]);
    cbz_mult(&(*gram)[0][1], &(*u)[0], &conj0);
    cbz_mult(&prod, &(*u)[1], &conj1);
    cbz_mult_ibz(&prod, &prod, &(alg->p));
    cbz_add(&(*gram)[0][1], &(*gram)[0][1], &prod);

    cbz_mult(&(*gram)[1][1], &(*v)[0], &conj0);
    cbz_mult(&prod, &(*v)[1], &conj1);
    cbz_mult_ibz(&prod, &prod, &(alg->p));
    cbz_add(&(*gram)[1][1], &(*gram)[1][1], &prod);

    mpz_mul_2exp((*gram)[0][0][0], (*gram)[0][0][0], 1);
    mpz_mul_2exp((*gram)[0][0][1], (*gram)[0][0][1], 1);
    mpz_mul_2exp((*gram)[0][1][0], (*gram)[0][1][0], 1);
    mpz_mul_2exp((*gram)[0][1][1], (*gram)[0][1][1], 1);
    mpz_mul_2exp((*gram)[1][0][0], (*gram)[1][0][0], 1);
    mpz_mul_2exp((*gram)[1][0][1], (*gram)[1][0][1], 1);
    mpz_mul_2exp((*gram)[1][1][0], (*gram)[1][1][0], 1);
    mpz_mul_2exp((*gram)[1][1][1], (*gram)[1][1][1], 1);

    cbz_finalize(&conj0); cbz_finalize(&conj1); cbz_finalize(&prod);
}

void cbz_lagrange_2x2(cbz_mat_2x2_t *uni_mat, cbz_mat_2x2_t *uv, cbz_mat_2x2_t *gram){

    cbz_set_uint32(&(*uni_mat)[0][0], 1, 0); cbz_set_uint32(&(*uni_mat)[0][1], 0, 0); cbz_set_uint32(&(*uni_mat)[1][0], 0, 0); cbz_set_uint32(&(*uni_mat)[1][1], 1, 0);

    if (ibz_cmp(&(*gram)[0][0][0], &(*gram)[1][1][0]) > 0){
        inter_gauss_lagrange_4x2_swap(uni_mat, uv, gram);
    }
    inter_gauss_lagrange_4x2_reduce(uni_mat, uv, gram);

    while ((ibz_cmp(&(*gram)[0][0][0], &(*gram)[1][1][0]) > 0)) {
        inter_gauss_lagrange_4x2_swap(uni_mat, uv, gram);
        inter_gauss_lagrange_4x2_reduce(uni_mat, uv, gram);
    }
}

void cbz_mod(cbz_t *r, const cbz_t *a, const ibz_t *mod){
    ibz_mod(&(*r)[0], &(*a)[0], mod);
    ibz_mod(&(*r)[1], &(*a)[1], mod);
}

void cbz_linear_comb_mod(cbz_t *r, const cbz_t *f1, const cbz_t *a, const cbz_t *f2, const cbz_t *b, const ibz_t *mod) {
    cbz_t prod1, prod2;
    cbz_init(&prod1); cbz_init(&prod2);

    cbz_mult(&prod1, f1, a);
    cbz_mult(&prod2, f2, b);
    cbz_add(r, &prod1, &prod2);
    cbz_mod(r, r, mod);

    cbz_finalize(&prod1); cbz_finalize(&prod2);
}

void cbz_vec_2_ibz_mod(cbz_vec_2_t *r, const cbz_vec_2_t *a, const ibz_t *mod){
    cbz_mod(&(*r)[0], &(*a)[0], mod);
    cbz_mod(&(*r)[1], &(*a)[1], mod);
}

int cbz_is_zero(const cbz_t *a){
    return ibz_is_zero(&(*a)[0]) && ibz_is_zero(&(*a)[1]);
}

int cbz_vec_2_is_zero(const cbz_vec_2_t *a){
    return cbz_is_zero(&(*a)[0]) && cbz_is_zero(&(*a)[1]);
}

void cbz_mat_2x2_scalar_mul(cbz_mat_2x2_t *res, const ibz_t *scalar, const cbz_mat_2x2_t *in_mat){
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            cbz_mult_ibz(&(*res)[i][j], &(*in_mat)[i][j], scalar);
        }   
    }
}

void cbz_hnf_mod(cbz_mat_2x2_t *hnf, int generator_number, const cbz_vec_2_t *generators, const ibz_t *mod, const quat_alg_t *alg) {
    assert(generator_number > 1);
    assert(!ibz_is_zero(mod));
    cbz_t tmp;
    cbz_vec_2_t augmented_gen[generator_number + 2];
    cbz_mat_2x2_t gram, uni_mat, uv;
    cbz_init(&tmp);
    cbz_mat_2x2_init(&gram); cbz_mat_2x2_init(&uni_mat); cbz_mat_2x2_init(&uv);
    for (int i = 0; i < generator_number; i++) {
        cbz_vec_2_init(&augmented_gen[i]);
        cbz_vec_2_ibz_mod(&augmented_gen[i], &generators[i], mod);
    }
    for (int i = generator_number; i < generator_number + 2; i++){
        cbz_vec_2_init(&augmented_gen[i]);
        cbz_set_ibz(&augmented_gen[i][1 - (i - generator_number)], mod, &ibz_const_zero);
    }

    cbz_set_uint32(&uv[0][1], 0, 0);
    cbz_set_uint32(&uv[1][1], 0, 0);
    for (int i = 1; i >= 0; i--){
        for (int j = 2-i; j < generator_number + 2 - i; j++) {
            cbz_copy(&uv[0][0], &augmented_gen[1-i][i]);
            cbz_copy(&uv[1][0], &augmented_gen[j][i]);
            cbz_gram_2x2(&gram, &uv, alg);
            cbz_lagrange_2x2(&uni_mat, &uv, &gram);
            cbz_copy(&augmented_gen[1-i][i], &uv[1][0]);
            cbz_copy(&augmented_gen[j][i], &uv[0][0]);
            for (int k = 0; k < i; k++){
                cbz_linear_comb_mod(&tmp, &uni_mat[1][0], &augmented_gen[1-i][k], &uni_mat[1][1], &augmented_gen[j][k], mod);
                cbz_linear_comb_mod(&augmented_gen[j][k], &uni_mat[0][0], &augmented_gen[1-i][k], &uni_mat[0][1], &augmented_gen[j][k], mod);
                cbz_copy(&augmented_gen[1-i][k], &tmp);
            }
        }
    }
    
    cbz_rounded_div(&tmp, &augmented_gen[0][0], &augmented_gen[0][0], &augmented_gen[1][0]);

    for (int i = 2; i < generator_number + 2; i++){
        assert(cbz_vec_2_is_zero(&augmented_gen[i]));
    }
    
    for (int i = 0; i < 2; i++){
        cbz_vec_2_copy(&(*hnf)[1-i], &augmented_gen[i]);
    }

    cbz_finalize(&tmp);
    cbz_mat_2x2_finalize(&gram); cbz_mat_2x2_finalize(&uni_mat); cbz_mat_2x2_finalize(&uv);
    for (int i = 0; i < generator_number + 2; i++){
        cbz_vec_2_finalize(&augmented_gen[i]);
    }
}

void
cbz_O0_lideal_create(quat_left_ideal_t *lideal,
                   const quat_alg_elem_t *x,
                   const ibz_t *N,
                   const quat_lattice_t *order,
                   const quat_alg_t *alg) {
    
    ibz_t mod;
    ibz_mat_4x4_t ON, Ox;
    cbz_vec_2_t generators[4];
    cbz_mat_2x2_t hnf, tmp;
    ibz_init(&mod);
    ibz_mat_4x4_init(&ON); ibz_mat_4x4_init(&Ox);
    cbz_mat_2x2_init(&hnf); cbz_mat_2x2_init(&tmp);
    for (int i = 0; i < 4; i++){
        cbz_vec_2_init(&generators[i]);
    }
    
    quat_lattice_mat_alg_coord_mul_without_hnf(&Ox, &(order->basis), &(x->coord), alg);
    ibz_mul(&(lideal->lattice.denom), &(order->denom), &(x->denom));
    ibz_copy(&(lideal->norm), N);
    lideal->parent_order = order;
    ibz_mat_4x4_scalar_mul(&ON, N, &order->basis);
    ibz_mat_4x4_scalar_mul(&ON, &(x->denom), &ON);
    
    cbz_mat_2x2_from_hnf_ibz_lat(&tmp, &Ox);
    cbz_vec_2_copy(&generators[0], &tmp[0]);
    cbz_vec_2_copy(&generators[1], &tmp[1]);
    cbz_mat_2x2_from_hnf_ibz_lat(&tmp, &ON);
    cbz_vec_2_copy(&generators[2], &tmp[0]);
    cbz_vec_2_copy(&generators[3], &tmp[1]);

    ibz_mul(&mod, N, &(lideal->lattice.denom));
    cbz_hnf_mod(&hnf, 4, generators, &mod, alg);
    ibz_lat_from_cbz_mat_2x2(&(lideal->lattice.basis), &hnf);
    quat_lattice_reduce_denom(&(lideal->lattice), &(lideal->lattice));

    assert(quat_lideal_norm_verify(lideal));
    
    ibz_finalize(&mod); 
    cbz_mat_2x2_finalize(&hnf); cbz_mat_2x2_finalize(&tmp);
    for (int i = 0; i < 4; i++){
        cbz_vec_2_finalize(&generators[i]);
    }
    ibz_mat_4x4_finalize(&ON); ibz_mat_4x4_finalize(&Ox);
}

void
cbz_mat_2x2_inv_given_denom(cbz_mat_2x2_t *inv, const ibz_t *denom, const cbz_mat_2x2_t *mat) {
    cbz_t const_denom, tmp, prod;
    cbz_mat_2x2_t work;
    cbz_init(&const_denom); cbz_init(&tmp); cbz_init(&prod);
    cbz_mat_2x2_init(&work);

    cbz_set_ibz(&const_denom, denom, &ibz_const_zero);
    cbz_div_floor_norm(&work[0][1], &tmp, &const_denom, &(*mat)[1][0]);
    assert(cbz_is_zero(&tmp));

    cbz_div_floor_norm(&work[1][0], &tmp, &const_denom, &(*mat)[0][1]);
    assert(cbz_is_zero(&tmp));

    cbz_mult(&prod, &(*mat)[0][0], &work[1][0]);
    cbz_neg(&prod, &prod);
    cbz_div_floor_norm(&work[1][1], &tmp, &prod, &(*mat)[1][0]);
    assert(cbz_is_zero(&tmp));

    for (int i = 0; i < 2; i++){
        cbz_vec_2_copy(&(*inv)[i], &work[i]);
    }

    cbz_finalize(&const_denom); cbz_finalize(&tmp); cbz_finalize(&prod);
    cbz_mat_2x2_finalize(&work);
}

void cbz_mat_2x2_transpose(cbz_mat_2x2_t *t_mat, const cbz_mat_2x2_t *mat) {
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            cbz_copy(&(*t_mat)[i][j], &(*mat)[i][j]);
        }
    }
    cbz_swap(&(*t_mat)[0][1], &(*t_mat)[1][0]);
}

void cbz_mat_2x2_conj_transpose(cbz_mat_2x2_t *t_mat, const cbz_mat_2x2_t *mat) {
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            cbz_conjugate(&(*t_mat)[i][j], &(*mat)[i][j]);
        }
    }
    cbz_swap(&(*t_mat)[0][1], &(*t_mat)[1][0]);
}

void cbz_mat_2x2_swap_columns(cbz_mat_2x2_t *t_mat, const cbz_mat_2x2_t *mat) {
    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            cbz_copy(&(*t_mat)[i][j], &(*mat)[i][j]);
        }
    }
    cbz_swap(&(*t_mat)[0][0], &(*t_mat)[0][1]);
    cbz_swap(&(*t_mat)[1][1], &(*t_mat)[1][0]);
}

void cbz_printf(const cbz_t *x){
    ibz_printf("%Zd, %Zd \n", (*x)[0], (*x)[1]);
}

void cbz_vec_2_printf(const cbz_vec_2_t *x){
    ibz_printf("[%Zd + %Zd*ii, %Zd + %Zd*ii] \n", (*x)[0][0], (*x)[0][1], (*x)[1][0], (*x)[1][1]);
}