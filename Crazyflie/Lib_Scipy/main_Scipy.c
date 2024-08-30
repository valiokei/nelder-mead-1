#include <stdio.h>
#include <math.h>

#define N 3 // Numero di variabili

double rosenbrock(double x[3])
{
    double a = 1.0, b = 100.0;
    return pow(a - x[0], 2) + b * pow(x[1] - x[0] * x[0], 2);
}

double rastrigin(double x[3])
{
    int n = 3;
    double A = 10.0;
    double sum = A * n;
    for (int i = 0; i < n; i++)
    {
        sum += x[i] * x[i] - A * cos(2 * M_PI * x[i]);
    }
    return sum;
}

double ackley(double x[3])
{
    double sum1 = 0.0;
    double sum2 = 0.0;
    for (int i = 0; i < 3; i++)
    {
        sum1 += x[i] * x[i];
        sum2 += cos(2 * M_PI * x[i]);
    }
    return -20.0 * exp(-0.2 * sqrt(sum1 / 3)) - exp(sum2 / 3) + 20 + M_E;
}

double sphere(double x[3])
{
    double sum = 0.0;
    for (int i = 0; i < 3; i++)
    {
        sum += x[i] * x[i];
    }
    return sum;
}

double beale(double x[3])
{
    double x0 = x[0];
    double x1 = x[1];
    return pow(1.5 - x0 + x0 * x1, 2) + pow(2.25 - x0 + x0 * x1 * x1, 2) + pow(2.625 - x0 + x0 * x1 * x1 * x1, 2);
}

static int function_type = 1; // 1: rosenbrock, 2: rastrigin, 3: ackley, 4: sphere, 5: beale

double objective_function(double x[3])
{
    // Esempio di funzione obiettivo
    // Definizione della variabile esterna

    // Selezione della funzione obiettivo in base alla variabile esterna
    double result;
    switch (function_type)
    {
    case 1:
        result = rosenbrock(x);
        break;
    case 2:
        result = rastrigin(x);
        break;
    case 3:
        result = ackley(x);
        break;
    case 4:
        result = sphere(x);
        break;
    case 5:
        result = beale(x);
        break;
    default:
        result = 0.0; // Funzione di default
        break;
    }

    return result;
}

void nelder_mead(double start[N], double result[N], int max_iter, double tol)
{
    double rho = 1.0, chi = 2.0, psi = 0.5, sigma = 0.5;
    double simplex[N + 1][N], fsimplex[N + 1];
    double centroid[N], xr[N], xe[N], xc[N];
    int i, j, iter = 0;
    double temp[N];

    // Inizializzazione del simplex
    for (i = 0; i < N + 1; ++i)
    {
        for (j = 0; j < N; ++j)
        {
            if (i == 0)
            {
                simplex[i][j] = start[j];
            }
            else
            {
                simplex[i][j] = start[j] + ((i == j + 1) ? 0.05 * start[j] : 0.00025);
            }
        }
        fsimplex[i] = objective_function(simplex[i]);
    }

    while (iter < max_iter)
    {
        // Ordinamento del simplex
        for (i = 0; i < N + 1; ++i)
        {
            for (j = i + 1; j < N + 1; ++j)
            {
                if (fsimplex[i] > fsimplex[j])
                {
                    double ftemp = fsimplex[i];
                    fsimplex[i] = fsimplex[j];
                    fsimplex[j] = ftemp;
                    for (int k = 0; k < N; ++k)
                    {
                        temp[k] = simplex[i][k];
                        simplex[i][k] = simplex[j][k];
                        simplex[j][k] = temp[k];
                    }
                }
            }
        }

        // Calcolo del centroide
        for (j = 0; j < N; ++j)
        {
            centroid[j] = 0.0;
            for (i = 0; i < N; ++i)
            {
                centroid[j] += simplex[i][j];
            }
            centroid[j] /= N;
        }

        // Operazione di riflessione
        for (j = 0; j < N; ++j)
        {
            xr[j] = centroid[j] + rho * (centroid[j] - simplex[N][j]);
        }
        double fxr = objective_function(xr);

        if (fxr < fsimplex[0])
        {
            // Operazione di espansione
            for (j = 0; j < N; ++j)
            {
                xe[j] = centroid[j] + chi * (xr[j] - centroid[j]);
            }
            double fxe = objective_function(xe);

            if (fxe < fxr)
            {
                for (j = 0; j < N; ++j)
                {
                    simplex[N][j] = xe[j];
                }
                fsimplex[N] = fxe;
            }
            else
            {
                for (j = 0; j < N; ++j)
                {
                    simplex[N][j] = xr[j];
                }
                fsimplex[N] = fxr;
            }
        }
        else if (fxr < fsimplex[N - 1])
        {
            for (j = 0; j < N; ++j)
            {
                simplex[N][j] = xr[j];
            }
            fsimplex[N] = fxr;
        }
        else
        {
            if (fxr < fsimplex[N])
            {
                for (j = 0; j < N; ++j)
                {
                    xc[j] = centroid[j] + psi * (xr[j] - centroid[j]);
                }
                double fxc = objective_function(xc);

                if (fxc <= fxr)
                {
                    for (j = 0; j < N; ++j)
                    {
                        simplex[N][j] = xc[j];
                    }
                    fsimplex[N] = fxc;
                }
                else
                {
                    for (i = 1; i < N + 1; ++i)
                    {
                        for (j = 0; j < N; ++j)
                        {
                            simplex[i][j] = simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j]);
                        }
                        fsimplex[i] = objective_function(simplex[i]);
                    }
                }
            }
            else
            {
                for (j = 0; j < N; ++j)
                {
                    xc[j] = centroid[j] - psi * (centroid[j] - simplex[N][j]);
                }
                double fxc = objective_function(xc);

                if (fxc < fsimplex[N])
                {
                    for (j = 0; j < N; ++j)
                    {
                        simplex[N][j] = xc[j];
                    }
                    fsimplex[N] = fxc;
                }
                else
                {
                    for (i = 1; i < N + 1; ++i)
                    {
                        for (j = 0; j < N; ++j)
                        {
                            simplex[i][j] = simplex[0][j] + sigma * (simplex[i][j] - simplex[0][j]);
                        }
                        fsimplex[i] = objective_function(simplex[i]);
                    }
                }
            }
        }

        iter++;

        // Verifica della convergenza
        double max_diff = 0.0;
        for (i = 0; i < N; ++i)
        {
            double diff = fabs(simplex[N][i] - simplex[0][i]);
            if (diff > max_diff)
            {
                max_diff = diff;
            }
        }

        if (max_diff < tol)
        {
            break;
        }
    }

    for (i = 0; i < N; ++i)
    {
        result[i] = simplex[0][i];
    }
}

int main()
{
    double start[N] = {1.10, 0.5, 4.0};
    double result[N];
    int max_iter = 200000;
    double tol = 1e-8;

    for (function_type = 1; function_type <= 5; function_type++)
    {
        nelder_mead(start, result, max_iter, tol);

        switch (function_type)
        {
        case 1:
            printf("Using Rosenbrock function\n");
            break;
        case 2:
            printf("Using Rastrigin function\n");
            break;
        case 3:
            printf("Using Ackley function\n");
            break;
        case 4:
            printf("Using Sphere function\n");
            break;
        case 5:
            printf("Using Beale function\n");
            break;
        }

        printf("Minimum found at: (%f, %f, %f)\n", result[0], result[1], result[2]);
        printf("\n");
    }

    return 0;
}