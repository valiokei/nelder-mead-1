#ifndef __NELDER_MEAD__
#define __NELDER_MEAD__

#include <stddef.h>

#define DIMENSIONS 3

// #include "model.h"

typedef float real;
typedef enum
{
    false,
    true
} bool;

/*
 * "Low-noise" squaring for arguments with no side-effects
 */
#define SQR(x) ((x) * (x))

/*
 * Function definition
 */
typedef struct Model model;

/*
 * Point definition
 *   x: n-dimensional array with point coordinates
 *   y: value of a function f applied to the coordinates x, y = f(x)
 */
typedef struct Point
{
    real x[DIMENSIONS];
    real y;
} point;

/*
 * Initialize function
 */
model *init_model(int idx, model *mdl);

/*
 * Return expected number of dimensions
 */
int dimensions(void);

/*
 * Cost function implementation
 *   model: model to optimize
 *   point: point where to evaluate the function
 */
void cost(const model *, point *);

#define BASE 10

/*
 * ANSI color codes
 * https://gist.github.com/RabaDabaDoba/145049536f815903c79944599c6f952a
 */
#define GRY "\x1B[1;30m"
#define RED "\x1B[1;31m"
#define GRN "\x1B[1;32m"
#define YLW "\x1B[1;33m"
#define BLU "\x1B[1;34m"
#define MGT "\x1B[1;35m"
#define CYN "\x1B[1;36m"
#define WHT "\x1B[1;37m"
#define NRM "\x1B[0m"
#define GRYBG "\x1B[0;100m"
#define REDBG "\x1B[0;101m"
#define GRNBG "\x1B[0;102m"
#define YLWBG "\x1B[0;103m"
#define BLUBG "\x1B[0;104m"
#define MGTBG "\x1B[0;105m"
#define CYNBG "\x1B[0;106m"
#define WHTBG "\x1B[0;107m"

/*
 * Optimizer settings
 */
typedef struct Optimset
{
    int precision;         // significant figures in floats/exponentials
    int format;            // fixed or exponential floating point format
    int verbose;           // toggle verbose output during minimization
    real tol_x;            // tolerance on the simplex solutions coordinates
    real tol_y;            // tolerance on the function value
    unsigned int max_iter; // maximum number of allowed iterations
    unsigned int max_eval; // maximum number of allowed function evaluations
    int adaptive;          // simplex updates reduced for dimension > 2
    real scale;            // scaling factor of initial simplex
} optimset;

/*
 * The "simplex" containing an array of n + 1 points each of dimension n
 */
typedef struct Simplex
{
    int n;
    unsigned int num_iter, num_eval;
    point vertices[DIMENSIONS + 1];
    point reflected;
    point expanded;
    point contracted;
    point centroid;
} simplex;

/*
 * "Simplex" or "Amoeba" optimizer
 */
void nelder_mead(const model *, const optimset *, simplex *, point *);

/*
 * Utility functions
 */
real distance(int, const point *, const point *);

int compare(const void *, const void *);

void sort(simplex *);

void init_simplex(int, real, const point *, simplex *smpl);

void update_centroid(simplex *);

void update_simplex(real, const point *, const point *, const point *,
                    const model *, simplex *, point *);

real tolerance_x(const simplex *);

real tolerance_y(const simplex *);

int terminated(const simplex *, const optimset *);

// point *init_point(int);

void copy_point(int, const point *, point *);

void print_point(int, const point *, int, int);

void print_value(real, int, int);

void print_verbose(int, const char *, ...);

void free_simplex(simplex *);

void free_point(point *);

#endif // __NELDER_MEAD__