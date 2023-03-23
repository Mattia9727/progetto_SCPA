#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

__global__ void matrix_vector_product(float *d_A, float *d_x, float *d_y, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        float sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += d_A[i * n + j] * d_x[j];
        }
        d_y[i] = sum;
    }
}

void matrix_vector_product_omp(float *A, float *x, float *y, int n, int threads) {
#pragma omp parallel for num_threads(threads)
    for (int i = 0; i < n; i++) {
        float sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i * n + j] * x[j];
        }
        y[i] = sum;
    }
}

int main() {
    int n = 10000; // dimensione matrice/vettore
    int threads = 4; // numero di thread per OpenMP

    // allocazione memoria sulla CPU
    float *A = (float*) malloc(n * n * sizeof(float));
    float *x = (float*) malloc(n * sizeof(float));
    float *y_cpu = (float*) malloc(n * sizeof(float));
    float *y_gpu = (float*) malloc(n * sizeof(float));

    // inizializzazione matrice/vettore
    for (int i = 0; i < n; i++) {
        x[i] = i;
        for (int j = 0; j < n; j++) {
            A[i * n + j] = i + j;
        }
    }

    // esecuzione prodotto matrice vettore su CPU con OpenMP
    double start_time = omp_get_wtime();
    matrix_vector_product_omp(A, x, y_cpu, n, threads);
    double end_time = omp_get_wtime();
    printf("Execution time on CPU with OpenMP: %f seconds\n
}