#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 1000

int main()
{
    int i, j;
    double A[N][N], x[N], y[N];

    // Inizializzazione delle matrici
    for (i = 0; i < N; i++) {
        x[i] = 1.0;
        y[i] = 0.0;
        for (j = 0; j < N; j++) {
            A[i][j] = rand() / (double)RAND_MAX;
        }
    }

    // Calcolo del prodotto matrice-vettore con OpenMP
#pragma omp parallel for
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            y[i] += A[i][j] * x[j];
        }
    }

    // Stampa del risultato
    for (i = 0; i < N; i++) {
        printf("%f ", y[i]);
    }
    printf("\n");

    return 0;
}
