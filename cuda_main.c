/*#include <stdio.h>
#include <stdlib.h>

#define N 1000
#define BLOCK_SIZE 256

__global__ void matrixVectorProduct(double *A, double *x, double *y)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        double sum = 0.0;
        for (int j = 0; j < N; j++) {
            sum += A[i*N+j] * x[j];
        }
        y[i] = sum;
    }
}

int main()
{
    int i, j;
    double *A, *x, *y; // matrice, vettori di input e di output
    double *d_A, *d_x, *d_y; // matrice, vettori di input e di output sulla GPU

    // Allocazione della memoria sulla CPU
    A = (double*)malloc(N*N*sizeof(double));
    x = (double*)malloc(N*sizeof(double));
    y = (double*)malloc(N*sizeof(double));

    // Inizializzazione delle matrici
    for (i = 0; i < N; i++) {
        x[i] = 1.0;
        y[i] = 0.0;
        for (j = 0; j < N; j++) {
            A[i*N+j] = rand() / (double)RAND_MAX;
        }
    }

    // Allocazione della memoria sulla GPU
    cudaMalloc(&d_A, N*N*sizeof(double));
    cudaMalloc(&d_x, N*sizeof(double));
    cudaMalloc(&d_y, N*sizeof(double));

    // Copia dei dati dalla CPU alla GPU
    cudaMemcpy(d_A, A, N*N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_x, x, N*sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, y, N*sizeof(double), cudaMemcpyHostToDevice);

    // Calcolo del prodotto matrice-vettore sulla GPU
    int numBlocks = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    matrixVectorProduct<<<numBlocks, BLOCK_SIZE>>>(d_A, d_x, d_y);

    // Copia dei dati dalla GPU alla CPU
    cudaMemcpy(y, d_y, N*sizeof(double), cudaMemcpyDeviceToHost);

    // Stampa del risultato
    for (i = 0; i < N; i++) {
        printf("%f ", y[i]);
    }
    printf("\n");

    // Liberazione della memoria
    free(A);
    free(x);
    free(y);
    cudaFree(d_A);
    cudaFree(d_x);
    cudaFree(d_y);

    return 0;
}
*/
