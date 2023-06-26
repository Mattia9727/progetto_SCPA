#ifndef _PRODUCTELLPACKGPUH_
#define _PRODUCTELLPACKGPUH_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/matrices/format/headers/ellpack.h"
#include "../src/matrices/headers/matrix_generator.h"
#include <cuda_runtime.h>  // For CUDA runtime API


#define NUM_THREADS 40

double* move_AS_on_gpu(ellpack_matrix* mat){
    double* AS;

    cudaError_t cudaStatus;

    cudaStatus = cudaMalloc((void**)&AS, sizeof(double) * mat->maxnz);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc AS failed! %s\n",cudaGetErrorString(cudaStatus));
    }
    cudaStatus = cudaMemcpy(AS, mat->AS,  sizeof(double) * mat->maxnz, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy AS failed! %s\n",cudaGetErrorString(cudaStatus));
    }

    return AS;
}

int* move_JA_on_gpu(ellpack_matrix* mat){
    int* JA;
    
    cudaError_t cudaStatus;

    cudaStatus = cudaMalloc((void**)&JA, sizeof(int) *mat->maxnz);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc JA failed! %s\n",cudaGetErrorString(cudaStatus));
    }
    cudaStatus = cudaMemcpy(JA, mat->JA,  sizeof(int) * mat->maxnz, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy JA failed! %s",cudaGetErrorString(cudaStatus));
    }

    return JA;
}


double* move_matrix_coeff_on_gpu(matrix* vector){
    double* coeff;

    cudaError_t cudaStatus;

    cudaStatus = cudaMalloc((void**)&coeff, sizeof(double) * vector->n * vector->m);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc coeff failed!\n");
    }
    cudaMemcpy(coeff, vector->coeff,  sizeof(double) * vector->n * vector->m, cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy coeff failed!\n");
    }
    return coeff;
}

void free_on_gpu(void* h_data){
    cudaError_t error;
    error = cudaFree(h_data);
    if( error != cudaSuccess ){
        fprintf(stderr, "cudaFree1 failed!\n");
    }
}

__global__ void optimized_cuda_ellpack_product(int m, int n, int maxnz, double* AS, int* JA, double* coeff, double* myRes){

    int idx = blockDim.x * blockIdx.x + threadIdx.x;

    if (idx < m*n) {                                    //Se l'id del thread supera il numero di elementi da calcolare ignoro il thread
        
        int start = (idx/n) * maxnz;                    //Per selezionare la riga di AS e JA
        int multiv_col = idx % n;                       //Per selezionare la colonna del multivettore
        double sum = 0;
        for (int i = 0; i < maxnz; i++) {
            int col = JA[start + i];
            double value = AS[start + i];
            sum += value * coeff[col+m*multiv_col];
        }
        
        myRes[idx] = sum;
    }

}

void print_gpu_matrix_d(double* matrix, int m, int n){
    double* app = (double*)malloc(sizeof(double)*m*n);
    cudaMemcpy(app,matrix,sizeof(double)*m*n,cudaMemcpyDeviceToHost);
    printf("MATRICE IN GPU\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", app[i * n + j]);
        }
        printf("\n");
    }
    free(app);
}

void print_gpu_matrix_i(int* matrix, int m, int n){
    int* app = (int*)malloc(sizeof(int)*m*n);
    cudaMemcpy(app,matrix,sizeof(int)*m*n,cudaMemcpyDeviceToHost);
    printf("MATRICE IN GPU\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", app[i * n + j]);
        }
        printf("\n");
    }
    free(app);
}

double optimized_cuda_ellpack_product_in(ellpack_matrix host_mat, matrix vector, matrix* result){

    printf("maxnz=%d\n", host_mat.maxnz);

    double* AS = move_AS_on_gpu(&host_mat);
    int* JA = move_JA_on_gpu(&host_mat);
    double* coeff = move_matrix_coeff_on_gpu(&vector);

    //print_ellpack_matrix(host_mat);
    //print_gpu_matrix_d(coeff,vector.m,vector.n);

    double* cuda_result;
    cudaMalloc((void**)&cuda_result, sizeof(double) * result->m * result->n);
    double* myRes = (double *)malloc(sizeof(double)*result->m*result->n);

    int blockX = int(((1024/vector.m)/32)*32), blockY = vector.m;  //gestione warp. Ogni blocco gestisce un numero di righe pari a blockX
    if (blockX==0) blockX++;
    int gridX = int(host_mat.m/blockX)+1;

    //printf("blockx=%d, blocky=%d, gridx=%d\n", blockX, blockY, gridX);

    dim3 gridSize(gridX);
    dim3 blockSize(blockX*blockY);    
    
    //dim3 gridSize(4);               //NUMERO DI BLOCCHI IN UNA GRID
    //dim3 blockSize(256);            //NUMERO DI THREAD IN UN BLOCCO

    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);// 0 - the default stream
    //optimized_cuda_ellpack_product(result->m, result->n, host_mat.maxnz, AS, JA, coeff, cuda_result);
    optimized_cuda_ellpack_product<<<gridSize,blockSize>>>(result->m, result->n, host_mat.maxnz, AS, JA, coeff, cuda_result);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);// 0 - the default stream
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("%f\n",time);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    result->coeff = myRes;
    free_on_gpu((void*)JA);
    free_on_gpu((void*)AS);
    free_on_gpu((void*)coeff);

    cudaMemcpy(result->coeff, cuda_result, sizeof(double) * result->n * result->m, cudaMemcpyDeviceToHost);
    
    free_on_gpu((void*)cuda_result);
    return (double)time;
}
#endif