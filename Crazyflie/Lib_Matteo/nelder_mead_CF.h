
#ifndef __NELDER_MEAD_CF__
#define __NELDER_MEAD_CF__

#include "Config.h"

#define DIMENSIONS 3

typedef long double real;

struct Model
{
    real V_Measured[NUM_ANCHORS];
};

/*
 * "Low-noise" squaring for arguments with no side-effects
 */
#define SQR(x) ((x) * (x))

/*
 * Function definition
 */
typedef struct Model model;
/*
 * Cost function implementation
 *   model: model to optimize
 *   point: point where to evaluate the function
 */

typedef struct Point
{
    real x[DIMENSIONS];
    real y;
} point;

typedef struct voltMeasurement_s
{

    real x[4];
    real y[4];
    real z[4];
    int anchorId[4];
    real measuredVolt[4];
    real resonanceFrequency[4];
    real stdDev[4];
    real GainValue;
    int Id_in_saturation;
} voltMeasurement_t;

/*
 * Initialize function
 */
void init_model(model mdl, int idxToTake);
void cost(const model *, point *);
/*
 * Point definition
 *   x: n-dimensional array with point coordinates
 *   y: value of a function f applied to the coordinates x, y = f(x)
 */

/*
 * Return expected number of dimensions
 */
int dimensions(void);

/*
 * Cost function implementation
 *   model: model to optimize
 *   point: point where to evaluate the function
 */

/*
 * Optimizer settings
 */
typedef struct Optimset
{
    int precision;         // significant figures in reals/exponentials
    int format;            // fixed or exponential realing point format
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

void free_simplex(simplex *);

void free_point(point *);

#endif // __NELDER_MEAD_CF__