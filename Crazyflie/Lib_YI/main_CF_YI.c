#include "nelder_mead_CF_YI.h"
#include <stdio.h>
#include <math.h>

// Funzione obiettivo
double objective_function(const Vec *vec)
{
    // Ad esempio, una semplice funzione quadratica: f(x, y) = x^2 + y^2
    return vec->val[0] * vec->val[0] + vec->val[1] * vec->val[1];
}

int main()
{
    // Imposta un semplice simplex iniziale per una funzione a 2 variabili
    unsigned int n = 2;
    Vec *simplex[3];
    double v0[] = {0.0, 0.0};
    double v1[] = {1.0, 0.0};
    double v2[] = {0.0, 1.0};

    simplex[0] = Vec_from_list(n, v0);
    simplex[1] = Vec_from_list(n, v1);
    simplex[2] = Vec_from_list(n, v2);

    NelderMead *nm = NelderMead_create(1e-8, 1000);
    Vec *result = NelderMead_minimize(nm, objective_function, simplex, n);

    printf("Minimum found at: (%e, %e)\n", result->val[0], result->val[1]);

    Vec_destroy(result);
    NelderMead_destroy(nm);

    return 0;
}