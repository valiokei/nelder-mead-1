
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "nelder_mead_ORIG.h"
#include "Config__ORIG.h"
#include "model_ORIGIN.h"

// parametri per la funzione di costo
// Definizione delle strutture e delle variabili globali
#define Num_TagPoses 10
// Parametri coil trasmittente (TX)
real Anchors[Number_of_Anchors][3] = {
    {-0.50f, -0.25f, 0.0f},
    {+0.50f, -0.25f, 0.0f},
    {+0.50f, +0.25f, 0.0f},
    {-0.50f, +0.25f, 0.0f}};
// Parametri coil mobile (TX)
const real h_tx = 0.0f; // altezza del nodo mobile in [mm]

real InputPoint[Num_TagPoses][3];

real Tag[Num_TagPoses][3];

const real f[Number_of_Anchors] = {213e3f, 203e3f, 193e3f, 183e3f};

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
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[0],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[1],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[2]};

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

real V_from_B(real *B_field, real *rx_versor, real resonanceFreq, real Gain)
{
    real dot_product_V = dot_product(B_field, rx_versor, 3);

    real V = Gain * fabsf(2.0f * PI * resonanceFreq * PI * RAY * RAY * N_WOUNDS * dot_product_V);
    return V;
}

void write_positions_to_file(float (*estimated_positions)[3], int num_positions)
{

    // Anchors File
    FILE *fileAnchors = fopen("Anchors_ORIG.csv", "w");
    if (fileAnchors == NULL)
    {
        perror("Impossibile aprire il fileAnchors");
        return;
    }

    // Scrive l'intestazione del fileAnchors CSV
    fprintf(fileAnchors, "AnchorX,AnchorY,AnchorZ\n");

    for (int anchorIdx = 0; anchorIdx < Number_of_Anchors; anchorIdx++)
    {
        fprintf(fileAnchors, "%f,%f,%f\n",
                (float)Anchors[anchorIdx][0], (float)Anchors[anchorIdx][1], (float)Anchors[anchorIdx][2]);
    }
    fclose(fileAnchors);

    // True Positions Tag File
    FILE *fileTrueTagPositions = fopen("TrueTagPositions_ORIG.csv", "w");
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
    FILE *fileEstimatedTagPositions = fopen("EstimatedTagPositions_ORIG.csv", "w");
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
                (float)estimated_positions[i][0], (float)estimated_positions[i][1], (float)estimated_positions[i][2]);
    }
    fclose(fileEstimatedTagPositions);

    // Input Points File
    FILE *fileInputPoint = fopen("InputPoint_ORIG.csv", "w");
    if (fileInputPoint == NULL)
    {
        perror("Impossibile aprire il fileEstimatedTagPositions");
        return;
    }
    fprintf(fileInputPoint, "InputPoint_x,InputPoint_y,InputPoint_z\n");
    for (int i = 0; i < num_positions; i++)
    {
        fprintf(fileInputPoint, "%f,%f,%f\n",
                (float)InputPoint[i][0], (float)InputPoint[i][1], (float)InputPoint[i][2]);
    }
    fclose(fileInputPoint);
}

optimset opt = {
    .precision = 100.0f,
    .format = 1,
    .verbose = 0,
    .tol_x = 1e-6f,
    .tol_y = 1e-6f,
    .max_iter = 5000.0f,
    .max_eval = 5000.0f,
    .adaptive = 0,
    .scale = 1.0e-5f};

struct Model
{
    real V_measured[Number_of_Anchors];
};
size_t dimensions()
{
    return PROBLEM_DIMENSION; // Tre componenti per la posizione (x, y, z)
}

model *init_model(int idxToTake)
{
    model *mdl = malloc(sizeof(model));
    ///============================= initialization of the stuff for the magnetic simulation ==============================

    // fill the model with the measured Volts

    real B_field[3];
    real V_measured[Number_of_Anchors];

    for (int anchorIdx = 0; anchorIdx < Number_of_Anchors; anchorIdx++)
    {
        get_B_field_for_a_Anchor(Anchors[anchorIdx], Tag[idxToTake], versore_spira_rx, B_field);
        V_measured[anchorIdx] = V_from_B(B_field, versore_spira_rx, f[anchorIdx], G_INA);
        mdl->V_measured[anchorIdx] = V_measured[anchorIdx];
    }

    return mdl;
}

void cost(const model *mdl, point *pnt)
{

    real costo = 0.0f;
    real V_predictedFromCurrentPnt[Number_of_Anchors];
    real B_field_vector[3];
    for (int anchorIdx = 0; anchorIdx < Number_of_Anchors; anchorIdx++)
    {

        get_B_field_for_a_Anchor(Anchors[anchorIdx], pnt->x, versore_spira_rx, B_field_vector);
        V_predictedFromCurrentPnt[anchorIdx] = V_from_B(B_field_vector, versore_spira_rx, f[anchorIdx], G_INA);
    }

    for (int i = 0; i < Number_of_Anchors; i++)
    {
        costo += powf(V_predictedFromCurrentPnt[i] - mdl->V_measured[i], 2);
    }

    pnt->y = costo;
}

int idxToTake = 0;

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
        const int n = PROBLEM_DIMENSION;
        // assert number of dimensions is correct
        // ASSERT(n == dimensions());

        // allocate input / output points
        // point *inp = init_point(n);
        // point *out = init_point(n);
        point *inp = init_point(n);
        point *out = init_point(n);

        // set input point coordinates from command args
        for (int i = 0; i < n; i++)
        {

            // True point
            // inp->x[i] = Tag[idxToTake][i];

            // point estimated at the previous iteration
            if (idxToTake == 0)
            {
                inp->x[i] = Tag[idxToTake][i];
            }
            else
            {
                inp->x[i] = out->x[i];
            }
            InputPoint[idxToTake][i] = inp->x[i];
        }

        // initialize model and simplex
        model *mdl = init_model(idxToTake + 1);

        simplex *smpl = init_simplex(n, opt.scale, inp);

        // optimize model mdl, using settings opt, starting
        // from simplex smpl, and save result to point out
        nelder_mead(mdl, &opt, smpl, out);
        cost(mdl, inp);
        printf("%s        Input%s ", WHT, NRM);
        print_point(n, inp, opt.precision, opt.format);
        printf("%s       Output%s ", WHT, NRM);

        // save the estimated position
        estimated_positions[idxToTake][0] = out->x[0];
        estimated_positions[idxToTake][1] = out->x[1];
        estimated_positions[idxToTake][2] = out->x[2];

        // compute the euclidean distance between the input and output points
        real euclidean_distance = sqrtf(pow(out->x[0] - Tag[idxToTake][0], 2) + pow(out->x[1] - Tag[idxToTake][1], 2) + pow(out->x[2] - Tag[idxToTake][2], 2));
        // printf("euclidean distance: %Lf\n", euclidean_distance);

        // print_point(n, inp, opt.precision, opt.format);

        // free memory and exit
        free_simplex(smpl);
        free_point(inp);
        free_point(out);
        free(mdl);
    }

    // write the files
    write_positions_to_file(estimated_positions, Num_TagPoses);

    // Liberazione della memoria allocata
    // free(x);
    // free(y);

    return 0;
}
