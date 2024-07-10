
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "nelder_mead.h"

// parametri per la funzione di costo
// Definizione delle strutture e delle variabili globali

#define Num_TagPoses 100
#define Num_Anchors 4

// Parametri coil trasmittente (TX)
real Anchors[Num_Anchors][3] = {
    {-0.50f, -0.25f, 0.0f},
    {+0.50f, -0.25f, 0.0f},
    {+0.50f, +0.25f, 0.0f},
    {-0.50f, +0.25f, 0.0f}};
// Parametri coil mobile (TX)
const real h_tx = 0.0f; // altezza del nodo mobile in [mm]

real InputPoint[Num_TagPoses][3];

real Tag[Num_TagPoses][3];
// real Tag[Num_TagPoses][3] = { // posizioni vere del mobile in [mm], poi convertite in [m]
//     {75.0f / 1000.0f, 0.00000001f / 1000.0f, h_tx},
//     {225.0f / 1000.0f, 0.00000001f / 1000.0f, h_tx},
//     {0.00000001f / 1000.0f, 75.0f / 1000.0f, h_tx},
//     {75.0f / 1000.0f, 75.0f / 1000.0f, h_tx},
//     {150.0f / 1000.0f, 75.0f / 1000.0f, h_tx},
//     {225.0f / 1000.0f, 75.0f / 1000.0f, h_tx},
//     {300.0f / 1000.0f, 75.0f / 1000.0f, h_tx},
//     {75.0f / 1000.0f, 150.0f / 1000.0f, h_tx},
//     {150.0f / 1000.0f, 150.0f / 1000.0f, h_tx},
//     {225.0f / 1000.0f, 150.0f / 1000.0f, h_tx},
//     {0.00000001f / 1000.0f, 225.0f / 1000.0f, h_tx},
//     {75.0f / 1000.0f, 225.0f / 1000.0f, h_tx},
//     {150.0f / 1000.0f, 225.0f / 1000.0f, h_tx},
//     {225.0f / 1000.0f, 225.0f / 1000.0f, h_tx},
//     {300.0f / 1000.0f, 225.0f / 1000.0f, h_tx},
//     {75.0f / 1000.0f, 300.0f / 1000.0f, h_tx},
//     {225.0f / 1000.0f, 300.0f / 1000.0f, h_tx}};

// Parametri coil riceventi (RX)
#define n_spire_rx 5                              // numero di spire in ciascun solenoide
#define raggio_spira_rx 0.019f                    // raggio di ciascun solenoide [m]
#define S_rx raggio_spira_rx *raggio_spira_rx *PI // area singola spira
#define G_INA 1000                                // guadagno INA
#define N_rx 16                                   // numero di coil RX

#define corrente_solenoide_tx 0.5f                  // Intensita di corrente che scorre in solenoide tx
#define raggio_spira_tx 0.019f                      // in [m]
#define S_tx raggio_spira_rx *raggio_spira_rx *M_PI // area singola spira
#define n_spire_tx 5
const real f[Num_Anchors] = {213e3f, 203e3f, 193e3f, 183e3f};

// #define RAY 0.019f
// #define N_WOUNDS 5.0f
// #define COIL_SURFACE (RAY * RAY * PI)
#define MU_0 1.25663706212e-06f
#define PI 3.14159265359f

real versore_spira_rx[3] = {0.0000000f, 0.0000000f, 1.0f};

real dot_product(real *a, real *b, int length)
{
    real result = 0.0f;
    for (int i = 0; i < length; i++)
    {
        result += a[i] * b[i];
    }
    return result;
}

real euclidean_distance(real *a, real *b, int length)
{
    real sum = 0.0f;
    for (int i = 0; i < length; i++)
    {
        real diff = a[i] - b[i];
        sum += diff * diff;
    }
    return sqrtf(sum);
}

void getversor(real *a, real *b, real *u, int length)
{
    real d = euclidean_distance(a, b, length);
    for (int i = 0; i < length; i++)
    {
        u[i] = (a[i] - b[i]) / d;
    }
}

void get_B_field_for_a_Anchor(real *anchor_pos,
                              real *tag_pos,
                              real *tag_or_versor,
                              real *B_field)
{

    real tx_rx_versor[3];
    getversor(anchor_pos, tag_pos, tx_rx_versor, 3);

    real magnetic_dipole_moment_tx[3] = {
        n_spire_rx * S_rx * corrente_solenoide_tx * tag_or_versor[0],
        n_spire_rx * S_rx * corrente_solenoide_tx * tag_or_versor[1],
        n_spire_rx * S_rx * corrente_solenoide_tx * tag_or_versor[2]};

    real tx_rx_distance = euclidean_distance(anchor_pos, tag_pos, 3);

    real dot_product_B_temp = dot_product(magnetic_dipole_moment_tx, tx_rx_versor, 3);

    real constant_Bfield_constant_1 = (MU_0 / (4.0f * PI)) / powf(tx_rx_distance, 3);
    real B_temp[3] = {
        constant_Bfield_constant_1 * (3.0f * dot_product_B_temp * tx_rx_versor[0] - magnetic_dipole_moment_tx[0]),
        constant_Bfield_constant_1 * (3.0f * dot_product_B_temp * tx_rx_versor[1] - magnetic_dipole_moment_tx[1]),
        constant_Bfield_constant_1 * (3.0f * dot_product_B_temp * tx_rx_versor[2] - magnetic_dipole_moment_tx[2])};

    // real magnetic_dipole_moment_tx_magnitude = euclidean_distance(magnetic_dipole_moment_tx, magnetic_dipole_moment_tx, 3);
    // compute the norm of the magnetic dipole moment
    // real magnetic_dipole_moment_tx_magnitude = sqrtf(magnetic_dipole_moment_tx[0] * magnetic_dipole_moment_tx[0] +
    //                                                   magnetic_dipole_moment_tx[1] * magnetic_dipole_moment_tx[1] +
    //                                                   magnetic_dipole_moment_tx[2] * magnetic_dipole_moment_tx[2]);

    // normalizing the B field (NOT NEEDED!!!!)
    B_field[0] = B_temp[0]; // * magnetic_dipole_moment_tx_magnitude;
    B_field[1] = B_temp[1]; // * magnetic_dipole_moment_tx_magnitude;
    B_field[2] = B_temp[2]; // * magnetic_dipole_moment_tx_magnitude;
}

real V_from_B(real *B_field, real *rx_versor, real resonanceFreq)
{
    real dot_product_V = dot_product(B_field, rx_versor, 3);

    real V = G_INA * fabsf(2.0f * PI * resonanceFreq * PI * raggio_spira_rx * raggio_spira_rx * n_spire_rx * dot_product_V);
    return V;
}

optimset opt = {
    .precision = 5.0f,
    .format = 0,
    .verbose = 0,
    .tol_x = 1e-3,
    .tol_y = 1e-3,
    .max_iter = 500,
    .max_eval = 500,
    .adaptive = 0,
    .scale = 1.0e-3f};

int problem_dimension = 3;

struct Model
{
    real V[Num_Anchors];
};

model *init_model(int idxToTake, model *mdl)
{
    // model *mdl = malloc(Num_Anchors * sizeof(model));
    ///============================= initialization of the stuff for the magnetic simulation ==============================

    // calculate the B field and the V for each coil in each positions
    // this is the ground truth
    real B_field[3];
    real V[Num_TagPoses][Num_Anchors];
    for (int posIdx = 0; posIdx < Num_TagPoses; posIdx++)
    {

        for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
        {
            get_B_field_for_a_Anchor(Anchors[anchorIdx], Tag[posIdx], versore_spira_rx, B_field);
            V[posIdx][anchorIdx] = V_from_B(B_field, versore_spira_rx, f[anchorIdx]);
        }
    }
    // print all V
    memcpy(mdl->V, V[idxToTake], Num_Anchors * sizeof(real));
    // printf("Tag[%d] = %Lf, %Lf, %Lf\n", idxToTake, Tag[idxToTake][0], Tag[idxToTake][1], Tag[idxToTake][2]);
    // for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
    // {
    //     printf("V[0][%d] = %Lf\n", anchorIdx, V[0][anchorIdx]);
    //     printf("mdl->V[%d] = %Lf\n", anchorIdx, mdl->V[anchorIdx]);
    // }

    // getchar();
    return mdl;
}

int dimensions()
{
    return 3; // Tre componenti per la posizione (x, y, z)
}

// Funzione di costo che funziona su un punto usando i V misurati da ogni ancora, questa funzione deve ottimizzare una posizione alla volta.
void cost(const model *mdl, point *pnt)
{
    // printf("pnt->x = %Lf, %Lf, %Lf\n", pnt->x[0], pnt->x[1], pnt->x[2]);
    // printf("Magnetic Cost Function\n");
    // printf("Input Point\n");
    // print_point(3, pnt, 9, 0);
    real costo = 0.0f;
    real V_model[Num_Anchors];
    real B_field_vector[3];
    for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
    {

        get_B_field_for_a_Anchor(Anchors[anchorIdx], pnt->x, versore_spira_rx, B_field_vector);
        // V_model[anchorIdx] = V_from_B(B_field_vector, versore_spira_rx, f0);

        // version with noise on the V_model
        V_model[anchorIdx] = V_from_B(B_field_vector, versore_spira_rx, f[anchorIdx]) + 0.001f * ((real)rand() / (real)RAND_MAX);
    }

    for (int i = 0; i < Num_Anchors; i++)
    {
        costo += powf(mdl->V[i] - V_model[i], 2);
    }

    pnt->y = costo;
}

int idxToTake = 0;

void write_positions_to_file(float (*estimated_positions)[3], int num_positions)
{

    // Anchors File
    FILE *fileAnchors = fopen("Anchors.csv", "w");
    if (fileAnchors == NULL)
    {
        perror("Impossibile aprire il fileAnchors");
        return;
    }

    // Scrive l'intestazione del fileAnchors CSV
    fprintf(fileAnchors, "AnchorX,AnchorY,AnchorZ\n");

    for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
    {
        fprintf(fileAnchors, "%f,%f,%f\n",
                (float)Anchors[anchorIdx][0], (float)Anchors[anchorIdx][1], (float)Anchors[anchorIdx][2]);
    }
    fclose(fileAnchors);

    // True Positions Tag File
    FILE *fileTrueTagPositions = fopen("TrueTagPositions.csv", "w");
    if (fileTrueTagPositions == NULL)
    {
        perror("Impossibile aprire il fileTrueTagPositions");
        return;
    }

    // Scrive l'intestazione del fileTrueTagPositions CSV
    fprintf(fileTrueTagPositions, "True_T_x,True_T_y,True_T_z\n");

    for (int posIdx = 0; posIdx < Num_TagPoses; posIdx++)
    {
        fprintf(fileTrueTagPositions, "%f,%f,%f\n",
                (float)Tag[posIdx][0], (float)Tag[posIdx][1], (float)Tag[posIdx][2]);
    }
    fclose(fileTrueTagPositions);

    // // Estimated Positions Tag File
    FILE *fileEstimatedTagPositions = fopen("EstimatedTagPositions.csv", "w");
    if (fileEstimatedTagPositions == NULL)
    {
        perror("Impossibile aprire il fileEstimatedTagPositions");
        return;
    }

    // Scrive l'intestazione del fileEstimatedTagPositions CSV
    fprintf(fileEstimatedTagPositions, "Estimated_T_x,Estimated_T_y,Estimated_T_z\n");

    for (int i = 0; i < num_positions; i++)
    {
        fprintf(fileEstimatedTagPositions, "%f,%f,%f\n",
                estimated_positions[i][0], estimated_positions[i][1], estimated_positions[i][2]);
    }
    fclose(fileEstimatedTagPositions);

    // Input Points File
    FILE *fileInputPoint = fopen("InputPoint.csv", "w");
    if (fileInputPoint == NULL)
    {
        perror("Impossibile aprire il fileEstimatedTagPositions");
        return;
    }
    fprintf(fileInputPoint, "InputPoint_x,InputPoint_y,InputPoint_z\n");
    for (int i = 0; i < num_positions; i++)
    {
        fprintf(fileInputPoint, "%f,%f,%f\n",
                InputPoint[i][0], InputPoint[i][1], InputPoint[i][2]);
    }
    fclose(fileInputPoint);
}

int main()
{

    real centro[2] = {0.0f, 0.0f}; // Centro della circonferenza
    real raggio = 0.3;             // Raggio della circonferenza
    int num_punti = Num_TagPoses;  // Numero di punti da generare
    // real *x, *y;
    real theta;
    int i;

    // Allocazione dinamica della memoria per le coordinate x e y
    // x = (real *)malloc(num_punti * sizeof(real));
    // y = (real *)malloc(num_punti * sizeof(real));
    real x[num_punti];
    real y[num_punti];

    // Generazione dei punti sulla circonferenza
    for (i = 0; i < num_punti; i++)
    {
        theta = 2 * PI * i / num_punti;         // Angolo per il punto corrente
        x[i] = centro[0] + raggio * cos(theta); // Coordinate x
        y[i] = centro[1] + raggio * sin(theta); // Coordinate y
    }
    for (i = 0; i < num_punti; i++)
    {
        Tag[i][0] = x[i];
        Tag[i][1] = y[i];
        Tag[i][2] = h_tx;
    }

    float estimated_positions[Num_TagPoses][3];

    // iterate over the tag positions
    for (int idxToTake = 0; idxToTake < Num_TagPoses; idxToTake++)
    {

        // ASSERT(opt.precision >= 3 && opt.precision <= 36);
        // ASSERT(opt.verbose == 0 || opt.verbose == 1);
        // ASSERT(opt.tol_x >= 1.0e-36f && opt.tol_x <= 1.0e-3f);
        // ASSERT(opt.tol_y >= 1.0e-36f && opt.tol_y <= 1.0e-3f);
        // ASSERT(opt.max_iter >= 1 && opt.max_iter <= 100000);
        // ASSERT(opt.max_eval >= 1 && opt.max_eval <= 100000);
        // ASSERT(opt.adaptive == 0 || opt.adaptive == 1);
        // ASSERT(opt.scale >= 1.0e-12f && opt.scale <= 1.0e3f);

        /// ======================== FUNZIONE DI COSTO ========================
        // questa funzione costo stima 1 volta, dovro iterare tutto per farla funzionare con putni multipli

        // read optimizer settings from command line and check values

        // infer number of dimension from command args
        const int n = problem_dimension;
        // assert number of dimensions is correct
        // ASSERT(n == dimensions());

        // allocate input / output points
        // point *inp = init_point(n);
        // point *out = init_point(n);
        point inp;
        point out;

        // set input point coordinates from command args
        for (int i = 0; i < n; i++)
        {
            // barycentre of the tag positions
            // real theta_0[3] = {0.15f, 0.15f, 0.2f};
            // inp.x[i] = 0.0;

            // gaussian noise to the input point applyied to the true point
            InputPoint[idxToTake][i] = Tag[idxToTake][i] + 0.1f * ((real)rand() / (real)RAND_MAX);
            inp.x[i] = InputPoint[idxToTake][i];

            // True point
            // inp->x[i] = Tag[idxToTake][i];
        }

        // initialize model and simplex
        model mdl;
        init_model(idxToTake, &mdl);

        simplex smpl;
        init_simplex(n, opt.scale, &inp, &smpl);

        // print information
        cost(&mdl, &inp);

        // optimize model mdl, using settings opt, starting
        // from simplex smpl, and save result to point out
        nelder_mead(&mdl, &opt, &smpl, &out);

        // printf("True point: [%Lf,%Lf,%Lf] \n", Tag[idxToTake][0], Tag[idxToTake][1], Tag[idxToTake][2]);
        // printf("starting point: [%Lf,%Lf,%Lf] \n", inp->x[0], inp->x[1], inp->x[2]);
        // printf("output point: [%Lf,%Lf,%Lf] \n", out->x[0], out->x[1], out->x[2]);

        // save the estimated position
        estimated_positions[idxToTake][0] = out.x[0];
        estimated_positions[idxToTake][1] = out.x[1];
        estimated_positions[idxToTake][2] = out.x[2];

        // compute the euclidean distance between the input and output points
        real euclidean_distance = sqrtf(pow(out.x[0] - Tag[idxToTake][0], 2) + pow(out.x[1] - Tag[idxToTake][1], 2) + pow(out.x[2] - Tag[idxToTake][2], 2));
        // printf("euclidean distance: %Lf\n", euclidean_distance);

        // print_point(n, inp, opt.precision, opt.format);

        // printf("%s        Input%s ", WHT, NRM);
        // print_point(n, inp, opt.precision, opt.format);
        // printf("%s       Output%s ", WHT, NRM);
        // print_point(n, out, opt.precision, opt.format);
        // printf("%s  Tolerance-x%s", WHT, NRM);
        // print_value(tolerance_x(smpl), opt.precision, opt.format);
        // printf("\n");
        // printf("%s  Tolerance-y%s", WHT, NRM);
        // print_value(tolerance_y(smpl), opt.precision, opt.format);
        // printf("\n");
        // printf("%s   Iterations%s %d\n", WHT, NRM, smpl->num_iter);
        // printf("%s  Evaluations%s %d\n", WHT, NRM, smpl->num_eval);

        // free memory and exit
        // free_simplex(&smpl);
        // free_point(&inp);
        // free_point(&out);
        // free(&mdl);
    }

    // write the files
    write_positions_to_file(estimated_positions, Num_TagPoses);

    // Liberazione della memoria allocata
    // free(x);
    // free(y);

    return 0;
}
