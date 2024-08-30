#ifndef __MAGNETIC_DECK_H__
#define __MAGNETIC_DECK_H__

#define PI 3.14159265359L

// -------------------------- ADC -------------------------------------------
// ADC DMA configuration
#define ARRAY_SIZE 2048

// 2^12 WHERE 12 IS THE NUMBER OF BITS OF THE MCU ADC = 4096
#define ADC_LEVELS 4096
#define ADC_MAX_VOLTAGE 3.0f
#define PCLK2 84e6f
#define ADC_PRESCALER 6.0f
// 12 from bit and 15 from the register value 12+15 = 27
#define ADC_Full_Sampling_Time 27.0f
#define Fc_ADC (PCLK2 / ADC_PRESCALER / ADC_Full_Sampling_Time)

#define ADC_Channel_Default ADC_Channel_3;
// -------------------------- DMA -------------------------------------------
#define DMA_IRQ DMA2_Stream4_IRQn
#define MY_DMA_Channel DMA_Channel_0
#define MY_DMA_Stream DMA2_Stream4

// IRQn_Type DMA_IRQ = DMA2_Stream4_IRQn;

// -------------------------- FFT -------------------------------------------
#define FFT_SIZE ARRAY_SIZE
#define BIN_SIZE (int)(Fc_ADC / FFT_SIZE)

// ------------------------ Measurement Model Params -------------------------------------------
#define Default_MagneticStandardDeviation 0.0001f
#define G_INA 700.0f

// -------------------------  adaptive std on Measured Voltage -------------------------------
#define UseAdaptiveSTD 0
#define window_size 25

// ------------------------ Calibration -------------------------------------------
#define CALIBRATION_TIC_VALUE 500.0f // Number of measurements to perform the calibration

// ------------------------System HZ-------------------------------------------
#define SYSTEM_HZ 100
#define SYSTEM_PERIOD_MS (1000 / SYSTEM_HZ)

// ------------------------ Anchors Parameters -------------------------------------------
#define NUM_ANCHORS 4

// Resonance Freqs Anchors in Hz
#define NeroResFreq 213e3
#define NeroIdx (int)(NeroResFreq / BIN_SIZE)
#define Nero_M -2.804
#define Nero_Q -2.635
#define Nero_Position_x -0.72f
#define Nero_Position_y -0.76f
#define Nero_Position_z +0.78f
#define Nero_Id 0

#define GialloResFreq 203e3
#define GialloIdx (int)(GialloResFreq / BIN_SIZE)
#define Giallo_M -2.887
#define Giallo_Q -2.629
#define Giallo_Position_x +0.64f
#define Giallo_Position_y +0.5f
#define Giallo_Position_z +0.78f
#define Giallo_Id 1

#define GrigioResFreq 193e3
#define GrigioIdx (int)(GrigioResFreq / BIN_SIZE)
#define Grigio_M -2.902
#define Grigio_Q -2.647
#define Grigio_Position_x -0.64f
#define Grigio_Position_y +0.5f
#define Grigio_Position_z +0.78f
#define Grigio_Id 2

#define RossoResFreq 183e3
#define RossoIdx (int)(RossoResFreq / BIN_SIZE)
#define Rosso_M -2.950
#define Rosso_Q -2.640
#define Rosso_Position_x +0.64f
#define Rosso_Position_y -0.5f
#define Rosso_Position_z +0.78f
#define Rosso_Id 3

// ------------------------ PHYSICAL COIL -------------------------------------------
#define RAY 0.019f
#define N_WOUNDS 5.0f
#define COIL_SURFACE (RAY * RAY * PI)
#define CURRENT 0.5f // This maybe can be improved
#define MU_0 1.25663706212e-06f

// ------------------------ Gain -------------------------------------------
// Potentiometer Params
#define R10 200.0f  // 200 Ohm
#define RW_2_7V 155 // ohm  155 OHM TYPICAL_DC_WIPER_RESISTANCE
#define RW_5_5V 100 // ohm

#define POTENTIOMETER_ADDR 0x2F // 0x94 WRITE 0x95 READ
#define POTENTIOMETER_BIT 7
// #define POTENTIOMETER_STEPS pow(2, POTENTIOMETER_BIT) // 7 bit --> 128
#define POTENTIOMETER_NUMBER_OF_STEPS (1 << POTENTIOMETER_BIT) - 1 // 127
#define POTENTIOMETER_FULL_SCALE_RAB 50E3                          // 50 kOhm
#define POTENTIOMETER_ANALOG_HW_DELAY_AFTER_SET 10
// in serie con l'INa c'e un operazionale che produce un guadagno di 10
#define DefaultPotentiometerValue G_INA / 10.0f
#define OpAmpGainValue 10.0f

// -------------------------- DAC -------------------------------------------
// DAC Reference Voltage Params
#define V_REF_CRAZYFLIE 3.0f
#define DAC_BIT 12
#define DAC_LEVELS (1 << DAC_BIT) // 4096
#define DAC_STEP (VREF / DAC_LEVELS)
#define DECK_DAC_I2C_ADDRESS 0x4C //
#define DAC_WRITE_LENGTH 2
#define DAC_ANALOG_HW_DELAY_AFTER_SET 10
#define V_DD 3.3f

// ========================== Function Definitions ==========================

// ----

// -------------------------- Measurement Model  -------------------------------------------

// void get_B_field_for_a_Anchor(float *anchor_pos,
//                               float *tag_pos,
//                               float *tag_or_versor,
//                               float *B_field);

// float V_from_B(float *B_field, float *rx_versor, float resonanceFreq, float Gain);

// // -------------------------- Math Utils -------------------------------------------

// float dot_product(float *a, float *b, int length);

// float euclidean_distance(float *a, float *b, int length);

// void getversor(float *a, float *b, float *u, int length);

// float computeSTD(float *data, int arrayDimension);

#endif // __MAGNETIC_DECK_H__