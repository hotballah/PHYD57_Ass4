/*
 * Filename:        Assn4.c
 * Description:     lifting-line theory of a full (test case) and a damaged wing of TU-154 airline.
 * Authors:          Abdilahi Mohamed, Ming Yan, Luke Shen, Yujun Liu
 * Last Updated:     2024-11-22
 *
 * Purpose:
 *   To predict the trajectory and orientation of the airliner after collision with a large tree
 *      and compare with accident reports. Confirm or debunk conspiracy theories surrounding this catastrophe.
 * 
 *
 *
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define N 1000
#define PI 3.14159265389793


const double wing_span = 10.0;
const double chord_length = 2.0;
const double U_inf = 15.0;
const double rho = 1.225;
const double angle_attack = 2;
const double zero_lift_angle = 0.0;


int main() {
    double y[N], Gamma[N], A[N][N], RHS[N], Lift[N];
    double dy = wing_span / (N - 1);
    double alpha = angle_attack * PI / 180.0;       
    double alpha_L0 = zero_lift_angle * PI / 180.0; 
    double dGamma, sum, L = 0.0;


    for (int i = 0; i < N; i++) {
        y[i] = -wing_span / 2 + i * dy; 
        Gamma[i] = 0.0;         
    }



for (int i = 0; i < N; i++) {
    RHS[i] = 2 * PI * U_inf * chord_length * (alpha - alpha_L0);
    for (int j = 0; j < N; j++) {
        if (i == j) {
            A[i][j] = 1.0; 
        } else {
            double delta_y = y[i] - y[j];
            A[i][j] = dy / (4.0 * PI * (delta_y * delta_y + 1)); // Biot-Savart induced velocity
        }
    }
}

for (int iter = 0; iter < 1000; iter++) {
    double max_error = 0.0;
    for (int i = 0; i < N; i++) {
        sum = 0.0;
        for (int j = 0; j < N; j++) {
            if (i != j) {
                sum += A[i][j] * Gamma[j];
            }
        }
        dGamma = (RHS[i] - sum) / A[i][i];
        max_error = fmax(max_error, fabs(dGamma - Gamma[i]));
        Gamma[i] = dGamma;
    }
    if (max_error < 1e-6) break; 
}


    
    for (int i = 0; i < N; i++) {
        Lift[i] = rho * U_inf * Gamma[i] * dy;
        L += Lift[i];
    }

    printf("Total Lift: %.2f N\n", L);

    for (int i = 0; i < N; i++){
        printf("%f,", Lift[i]);
    }
    return 0;
}
