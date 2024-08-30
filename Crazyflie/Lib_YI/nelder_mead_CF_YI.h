#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#ifdef __cplusplus
extern "C"
{
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

    typedef struct
    {
        unsigned int n;
        double *val;
    } Vec;

    Vec *Vec_create(unsigned int n);
    Vec *Vec_from_list(unsigned int n, double *values);
    Vec *Vec_copy(const Vec *src);
    double Vec_get(const Vec *vec, unsigned int idx);
    void Vec_set(Vec *vec, unsigned int idx, double value);
    void Vec_add(Vec *vec, const Vec *other);
    void Vec_sub(Vec *vec, const Vec *other);
    void Vec_div(Vec *vec, double scalar);
    void Vec_scale(Vec *vec, double scalar);
    unsigned int Vec_size(const Vec *vec);
    void Vec_destroy(Vec *vec);

    double inner_product(const Vec *a, const Vec *b);

    typedef struct
    {
        double tol;
        unsigned int max_iter;
        double alpha;
        double gamma;
        double rho;
        double sigma;
    } NelderMead;

    NelderMead *NelderMead_create(double tol, unsigned int max_iter);
    void NelderMead_destroy(NelderMead *nm);
    Vec *NelderMead_minimize(NelderMead *nm, double (*f)(const Vec *), Vec **simplex, unsigned int n);

#ifdef __cplusplus
}
#endif

#endif
