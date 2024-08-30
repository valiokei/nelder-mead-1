
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include "nelder_mead_CF.h"
#include "Config.h"
#include "time.h"

// parametri per la funzione di costo
// Definizione delle strutture e delle variabili globali

#define Num_TagPoses 10
#define Num_Anchors 4

// Parametri coil trasmittente (TX)
real Anchors[Num_Anchors][3] = {
    {-0.50f, -0.25f, 0.0f},
    {+0.50f, -0.25f, 0.0f},
    {+0.50f, +0.25f, 0.0f},
    {-0.50f, +0.25f, 0.0f}};
// Parametri coil mobile (TX)
const real h_tx = 0.0f; // altezza del nodo mobile in [mm]
// const real frequencies[NUM_ANCHORS] = {213e3f, 203e3f, 193e3f, 183e3f};
const real frequencies[NUM_ANCHORS] = {213e3f, 203e3f, 193e3f, 183e3f};

real InputPoint[Num_TagPoses][3];

real Tag[Num_TagPoses][3];

volatile real MeasuredVoltages_calibrated[4] = {0.0f, 0.0f, 0.0f, 0.0f};
static real estimated_position[3] = {0.0f, 0.0f, 0.0f};

real versore_spira_rx[3] = {0.0000f, 0.0f, 1.0f};

// --------------------------- Math Utils Functions ---------------------------
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
    return sqrtl(sum);
}

void getversor(real *a, real *b, real *u, int length)
{
    real d = euclidean_distance(a, b, length);
    for (int i = 0; i < length; i++)
    {
        u[i] = (a[i] - b[i]) / d;
    }
}

real computeSTD(real *data, int arrayDimension)
{
    real sum = 0.0;
    real mean = 0.0;
    real standardDeviation = 0.0;

    int i;

    for (i = 0; i < arrayDimension; i++)
    {
        sum += data[i];
    }

    mean = sum / arrayDimension;

    for (i = 0; i < arrayDimension; i++)
        standardDeviation += powl(data[i] - mean, 2);

    return sqrtf(standardDeviation / arrayDimension);
}

// ---------------- Measurement Model Functions ------------------------------
void get_B_field_for_a_Anchor(real *anchor_pos,
                              real *tag_pos,
                              real *tag_or_versor,
                              real *B_field)
{
    real tx_rx_versor[3];
    getversor(anchor_pos, tag_pos, tx_rx_versor, 3);
    // DEBUG_PRINT("tx_rx_versor = %f\n", tx_rx_versor);

    real magnetic_dipole_moment_tx[3] = {
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[0],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[1],
        N_WOUNDS * COIL_SURFACE * CURRENT * tag_or_versor[2]};
    // DEBUG_PRINT("magnetic_dipole_moment_tx = %f\n", magnetic_dipole_moment_tx);

    real tx_rx_distance = euclidean_distance(anchor_pos, tag_pos, 3);

    // DEBUG_PRINT("tx_rx_distance = %f\n", tx_rx_distance);

    real dot_product_B_temp = dot_product(magnetic_dipole_moment_tx, tx_rx_versor, 3);

    real constant_Bfield_constant_1 = (MU_0 / (4.0f * PI)) / powl(tx_rx_distance, 3);
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
static point inp;
static point out;
static model mdl;
static simplex smpl;

static real magnetiSTDNoise = 0;
const int problem_dimension = 3;
optimset opt = {

    // ASSERT(opt.precision >= 3 && opt.precision <= 36);
    // ASSERT(opt.verbose == 0 || opt.verbose == 1);
    // ASSERT(opt.tol_x >= 1.0e-36f && opt.tol_x <= 1.0e-3f);
    // ASSERT(opt.tol_y >= 1.0e-36f && opt.tol_y <= 1.0e-3f);
    // ASSERT(opt.max_iter >= 1 && opt.max_iter <= 100000);
    // ASSERT(opt.max_eval >= 1 && opt.max_eval <= 100000);
    // ASSERT(opt.adaptive == 0 || opt.adaptive == 1);
    // ASSERT(opt.scale >= 1.0e-12f && opt.scale <= 1.0e3f);

    .precision = 100.0f,
    .format = 0,
    .verbose = 0,
    .tol_x = 1e-6f,
    .tol_y = 1e-6f,
    .max_iter = 5000.0f,
    .max_eval = 5000.0f,
    .adaptive = 0,
    .scale = 1.0e-5f};

void init_model(model mdl, int idxToTake)
{
    ///============================= initialization of the stuff for the magnetic simulation ==============================

    // fill the model with the measured Volts

    real B_field[3];
    real V_measured[NUM_ANCHORS];

    for (int anchorIdx = 0; anchorIdx < NUM_ANCHORS; anchorIdx++)
    {
        get_B_field_for_a_Anchor(Anchors[anchorIdx], Tag[idxToTake], versore_spira_rx, B_field);
        V_measured[anchorIdx] = V_from_B(B_field, versore_spira_rx, frequencies[anchorIdx], G_INA);
        mdl.V_Measured[anchorIdx] = V_measured[anchorIdx];
    }
}

void cost(const model *mdl, point *pnt)
{

    real costo = 0.0f;
    real V_predictedFromCurrentPnt[NUM_ANCHORS];
    real B_field_vector[3];
    for (int anchorIdx = 0; anchorIdx < NUM_ANCHORS; anchorIdx++)
    {

        get_B_field_for_a_Anchor(Anchors[anchorIdx], pnt->x, versore_spira_rx, B_field_vector);
        V_predictedFromCurrentPnt[anchorIdx] = V_from_B(B_field_vector, versore_spira_rx, frequencies[anchorIdx], G_INA);
    }

    for (int i = 0; i < NUM_ANCHORS; i++)
    {
        costo += powl(V_predictedFromCurrentPnt[i] - mdl->V_Measured[i], 2);
    }

    pnt->y = costo;
}

int dimensions()
{
    return problem_dimension; // Tre componenti per la posizione (x, y, z)
}

void write_positions_to_file(real (*estimated_positions)[3], int num_positions)
{

    // Anchors File
    FILE *fileAnchors = fopen("Anchors_CF.csv", "w");
    if (fileAnchors == NULL)
    {
        perror("Impossibile aprire il fileAnchors");
        return;
    }

    // Scrive l'intestazione del fileAnchors CSV
    fprintf(fileAnchors, "AnchorX,AnchorY,AnchorZ\n");

    for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
    {
        fprintf(fileAnchors, "%Lf,%Lf,%Lf\n",
                (real)Anchors[anchorIdx][0], (real)Anchors[anchorIdx][1], (real)Anchors[anchorIdx][2]);
    }
    fclose(fileAnchors);

    // True Positions Tag File
    FILE *fileTrueTagPositions = fopen("TrueTagPositions_CF.csv", "w");
    if (fileTrueTagPositions == NULL)
    {
        perror("Impossibile aprire il fileTrueTagPositions");
        return;
    }

    // Scrive l'intestazione del fileTrueTagPositions CSV
    fprintf(fileTrueTagPositions, "True_T_x,True_T_y,True_T_z\n");

    for (int posIdx = 0; posIdx < Num_TagPoses; posIdx++)
    {
        fprintf(fileTrueTagPositions, "%Lf,%Lf,%Lf\n",
                (real)Tag[posIdx][0], (real)Tag[posIdx][1], (real)Tag[posIdx][2]);
    }
    fclose(fileTrueTagPositions);

    // // Estimated Positions Tag File
    FILE *fileEstimatedTagPositions = fopen("EstimatedTagPositions_CF.csv", "w");
    if (fileEstimatedTagPositions == NULL)
    {
        perror("Impossibile aprire il fileEstimatedTagPositions");
        return;
    }

    // Scrive l'intestazione del fileEstimatedTagPositions CSV
    fprintf(fileEstimatedTagPositions, "Estimated_T_x,Estimated_T_y,Estimated_T_z\n");

    for (int i = 0; i < num_positions; i++)
    {
        fprintf(fileEstimatedTagPositions, "%Lf,%Lf,%Lf\n",
                estimated_positions[i][0], estimated_positions[i][1], estimated_positions[i][2]);
    }
    fclose(fileEstimatedTagPositions);

    // Input Points File
    FILE *fileInputPoint = fopen("InputPoint_CF.csv", "w");
    if (fileInputPoint == NULL)
    {
        perror("Impossibile aprire il fileEstimatedTagPositions");
        return;
    }
    fprintf(fileInputPoint, "InputPoint_x,InputPoint_y,InputPoint_z\n");
    for (int i = 0; i < num_positions; i++)
    {
        fprintf(fileInputPoint, "%Lf,%Lf,%Lf\n",
                InputPoint[i][0], InputPoint[i][1], InputPoint[i][2]);
    }
    fclose(fileInputPoint);
}

int main()
{

    srand(time(NULL));

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

    // shift the first row of Tag to the last row
    real temp[3];
    temp[0] = Tag[0][0];
    temp[1] = Tag[0][1];
    temp[2] = Tag[0][2];

    for (int i = 0; i < Num_TagPoses - 1; i++)
    {
        Tag[i][0] = Tag[i + 1][0];
        Tag[i][1] = Tag[i + 1][1];
        Tag[i][2] = Tag[i + 1][2];
    }

    Tag[Num_TagPoses - 1][0] = temp[0];
    Tag[Num_TagPoses - 1][1] = temp[1];
    Tag[Num_TagPoses - 1][2] = temp[2];

    real estimated_positions[Num_TagPoses][3];

    // first point is the true point
    estimated_position[0] = Tag[0][0];
    estimated_position[1] = Tag[0][1];
    estimated_position[2] = Tag[0][2];
    // iterate over the tag positions
    for (int idxToTake = 0; idxToTake < Num_TagPoses; idxToTake++)
    {

        printf("TAG POSITION:  (%Lf,%Lf,%Lf)\n", Tag[idxToTake][0], Tag[idxToTake][1], Tag[idxToTake][2]);

        // filling  voltMeasurement_t *voltAnchor
        // static voltMeasurement_t voltAnchor;

        // // iterating over the anchors and fill voltAnchor
        // for (int anchorIdx = 0; anchorIdx < Num_Anchors; anchorIdx++)
        // {
        //     voltAnchor.anchorId[anchorIdx] = anchorIdx;
        //     voltAnchor.x[anchorIdx] = Anchors[anchorIdx][0];
        //     voltAnchor.y[anchorIdx] = Anchors[anchorIdx][1];
        //     voltAnchor.z[anchorIdx] = Anchors[anchorIdx][2];
        //     voltAnchor.resonanceFrequency[anchorIdx] = frequencies[anchorIdx];
        //     voltAnchor.GainValue = G_INA;

        //     // computing the Bfield for this anchor in the real position of the tag
        //     real B_field[3];
        //     get_B_field_for_a_Anchor(Anchors[anchorIdx], Tag[idxToTake + 1], versore_spira_rx, B_field);

        //     // computing the V for this anchor in the real position of the tag
        //     real MeasuredVolt = V_from_B(B_field, versore_spira_rx, voltAnchor.resonanceFrequency[anchorIdx], voltAnchor.GainValue);
        //     // adding gaussian noise to the measured voltage
        //     voltAnchor.measuredVolt[anchorIdx] = MeasuredVolt; // + magnetiSTDNoise * ((real)rand() / (real)RAND_MAX);
        // }

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

        // set input point coordinates from command args
        for (int i = 0; i < problem_dimension; i++)
        {
            // point estimated at the previous iteration
            if (idxToTake == 0)
            {
                inp.x[i] = Tag[idxToTake][i];
            }
            else
            {
                inp.x[i] = out.x[i];
            }
            InputPoint[idxToTake][i] = inp.x[i];
        }

        // initialize model and simplex
        real tag_or_versor[3] = {0.0, 0.0, 1.0};

        init_model(mdl, idxToTake + 1);

        init_simplex(n, opt.scale, &inp, &smpl);

        // optimize model mdl, using settings opt, starting
        // from simplex smpl, and save result to point out
        nelder_mead(&mdl, &opt, &smpl, &out);
        cost(&mdl, &inp);

        // save the estimated position
        // estimated_position[0] = out.x[0];
        // estimated_position[1] = out.x[1];
        // estimated_position[2] = out.x[2];
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
