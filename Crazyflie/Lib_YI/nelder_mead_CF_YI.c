#include "nelder_mead_CF_YI.h"

Vec *Vec_create(unsigned int n)
{
    Vec *vec = (Vec *)malloc(sizeof(Vec));
    vec->n = n;
    vec->val = (double *)calloc(n, sizeof(double));
    return vec;
}

Vec *Vec_from_list(unsigned int n, double *values)
{
    Vec *vec = (Vec *)malloc(sizeof(Vec));
    vec->n = n;
    vec->val = (double *)malloc(n * sizeof(double));
    memcpy(vec->val, values, n * sizeof(double));
    return vec;
}

Vec *Vec_copy(const Vec *src)
{
    return Vec_from_list(src->n, src->val);
}

double Vec_get(const Vec *vec, unsigned int idx)
{
    if (idx >= vec->n)
    {
        fprintf(stderr, "Element access out of range\n");
        exit(EXIT_FAILURE);
    }
    return vec->val[idx];
}

void Vec_set(Vec *vec, unsigned int idx, double value)
{
    if (idx >= vec->n)
    {
        fprintf(stderr, "Element access out of range\n");
        exit(EXIT_FAILURE);
    }
    vec->val[idx] = value;
}

void Vec_add(Vec *vec, const Vec *other)
{
    if (vec->n != other->n)
    {
        fprintf(stderr, "Vector lengths do not match\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < vec->n; i++)
    {
        vec->val[i] += other->val[i];
    }
}

void Vec_sub(Vec *vec, const Vec *other)
{
    if (vec->n != other->n)
    {
        fprintf(stderr, "Vector lengths do not match\n");
        exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < vec->n; i++)
    {
        vec->val[i] -= other->val[i];
    }
}

void Vec_div(Vec *vec, double scalar)
{
    for (unsigned int i = 0; i < vec->n; i++)
    {
        vec->val[i] /= scalar;
    }
}

void Vec_scale(Vec *vec, double scalar)
{
    for (unsigned int i = 0; i < vec->n; i++)
    {
        vec->val[i] *= scalar;
    }
}

unsigned int Vec_size(const Vec *vec)
{
    return vec->n;
}

void Vec_destroy(Vec *vec)
{
    free(vec->val);
    free(vec);
}

double inner_product(const Vec *a, const Vec *b)
{
    if (a->n != b->n)
    {
        fprintf(stderr, "Vector lengths do not match\n");
        exit(EXIT_FAILURE);
    }
    double sum = 0;
    for (unsigned int i = 0; i < a->n; i++)
    {
        sum += a->val[i] * b->val[i];
    }
    return sum;
}

NelderMead *NelderMead_create(double tol, unsigned int max_iter)
{
    NelderMead *nm = (NelderMead *)malloc(sizeof(NelderMead));
    nm->tol = tol;
    nm->max_iter = max_iter;
    nm->alpha = 1.0;
    nm->gamma = 2.0;
    nm->rho = 0.5;
    nm->sigma = 0.5;
    return nm;
}

void NelderMead_destroy(NelderMead *nm)
{
    free(nm);
}

double norm(Vec **simplex, unsigned int size)
{
    double sum = 0;
    for (unsigned int i = 1; i < size; i++)
    {
        Vec *diff = Vec_copy(simplex[i]);
        Vec_sub(diff, simplex[0]);
        sum += inner_product(diff, diff);
        Vec_destroy(diff);
    }
    return sqrt(sum / size);
}

Vec *NelderMead_minimize(NelderMead *nm, double (*f)(const Vec *), Vec **simplex, unsigned int n)
{
    Vec *centroid = Vec_create(n);
    for (unsigned int iter = 0; iter < nm->max_iter; iter++)
    {
        // Sort simplex vertices by function values
        for (unsigned int i = 0; i < n; i++)
        {
            for (unsigned int j = i + 1; j < n + 1; j++)
            {
                if (f(simplex[j]) < f(simplex[i]))
                {
                    Vec *temp = simplex[i];
                    simplex[i] = simplex[j];
                    simplex[j] = temp;
                }
            }
        }

        // Compute the centroid of the best n points
        for (unsigned int i = 0; i < n; i++)
        {
            Vec_add(centroid, simplex[i]);
        }
        Vec_div(centroid, n);

        // Reflection
        Vec *reflected = Vec_copy(centroid);
        Vec *diff = Vec_copy(centroid);
        Vec_sub(diff, simplex[n]);
        Vec_scale(diff, nm->alpha);
        Vec_add(reflected, diff);
        Vec_destroy(diff);

        if (f(reflected) < f(simplex[0]))
        {
            // Expansion
            Vec *expanded = Vec_copy(centroid);
            diff = Vec_copy(reflected);
            Vec_sub(diff, centroid);
            Vec_scale(diff, nm->gamma);
            Vec_add(expanded, diff);
            Vec_destroy(diff);

            if (f(expanded) < f(reflected))
            {
                Vec_destroy(simplex[n]);
                simplex[n] = expanded;
            }
            else
            {
                Vec_destroy(simplex[n]);
                simplex[n] = reflected;
                Vec_destroy(expanded);
            }
        }
        else if (f(reflected) < f(simplex[n - 1]))
        {
            Vec_destroy(simplex[n]);
            simplex[n] = reflected;
        }
        else
        {
            Vec *contracted;
            if (f(reflected) < f(simplex[n]))
            {
                // Outside contraction
                contracted = Vec_copy(centroid);
                diff = Vec_copy(reflected);
                Vec_sub(diff, centroid);
                Vec_scale(diff, nm->rho);
                Vec_add(contracted, diff);
                Vec_destroy(diff);
            }
            else
            {
                // Inside contraction
                contracted = Vec_copy(centroid);
                diff = Vec_copy(simplex[n]);
                Vec_sub(diff, centroid);
                Vec_scale(diff, nm->rho);
                Vec_add(contracted, diff);
                Vec_destroy(diff);
            }

            if (f(contracted) < f(simplex[n]))
            {
                Vec_destroy(simplex[n]);
                simplex[n] = contracted;
            }
            else
            {
                // Shrink
                for (unsigned int i = 1; i < n + 1; i++)
                {
                    diff = Vec_copy(simplex[i]);
                    Vec_sub(diff, simplex[0]);
                    Vec_scale(diff, nm->sigma);
                    Vec_add(simplex[i], diff);
                    Vec_destroy(diff);
                }
                Vec_destroy(contracted);
            }
            Vec_destroy(reflected);
        }

        if (norm(simplex, n + 1) < nm->tol)
        {
            break;
        }
    }

    // Sort simplex vertices by function values
    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = i + 1; j < n + 1; j++)
        {
            if (f(simplex[j]) < f(simplex[i]))
            {
                Vec *temp = simplex[i];
                simplex[i] = simplex[j];
                simplex[j] = temp;
            }
        }
    }

    Vec *result = Vec_copy(simplex[0]);
    for (unsigned int i = 0; i < n + 1; i++)
    {
        Vec_destroy(simplex[i]);
    }
    Vec_destroy(centroid);
    return result;
}
