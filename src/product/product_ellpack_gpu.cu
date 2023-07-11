#ifndef _PRODUCTELLPACKGPUH_
#define _PRODUCTELLPACKGPUH_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "../src/matrices/format/headers/ellpack.h"
#include "headers/product_ellpack.h"
#include "../src/matrices/headers/matrix_generator.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

#define MAX_NNZ_PER_WG 6144
#define MAX_BLOCK_THREADS 1024
#define MAX_GRID_SIZE 65536
#define WARP_SIZE 32

#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})


void print_gpu_matrix_d(double* matrix, int m, int n){
    double* app = (double*)malloc(sizeof(double)*m*n);
    checkCudaErrors(cudaMemcpy(app,matrix,sizeof(double)*m*n,cudaMemcpyDeviceToHost));
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
    checkCudaErrors(cudaMemcpy(app,matrix,sizeof(int)*m*n,cudaMemcpyDeviceToHost));
    printf("MATRICE IN GPU\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%d ", app[i * n + j]);
        }
        printf("\n");
    }
    free(app);
}


void print_gpu_matrix_l(long* matrix, long dim){
    long* app = (long*)malloc(sizeof(long)*dim);
    checkCudaErrors(cudaMemcpy(app,matrix,sizeof(long)*dim,cudaMemcpyDeviceToHost));
    printf("MATRICE IN GPU\n");
    for (long i = 0; i < dim; i++) {
        printf("%ld ", app[i]);
    }
    free(app);
}


__global__ void optimized_cuda_ellpack_product(int m, int n, int maxnz, double* AS, long* JA, double* coeff, double* myRes){

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


__global__ void optimized_cuda_h_ellpack_product(int m, int n, int* maxnz, double* AS, int* JA, int* hackOffsets, int hackSize, int numMatrix, int matDim, double* coeff, double* myRes){
    int bid = blockIdx.x;
    int tid = threadIdx.x;
    
    //Definisco vettori in shared memory che conterranno gli elementi delle "sotto-sottomatrici" di AS e JA
    __shared__ double vals[4096];
    __shared__ int cols[4096];

    int first = 0;
    int last = 0;
    if(bid == 0 && tid == 0) printf("%d\n",numMatrix);
    double res=0;
    for (int submatIdx = bid; submatIdx<numMatrix; submatIdx += gridDim.x){
        int counter = 0;
        for (first = hackOffsets[submatIdx]; first < hackOffsets[submatIdx + 1]; first += (int)((4096.0/maxnz[submatIdx]))*maxnz[submatIdx]){
            last = min(first + (int)((4096.0/maxnz[submatIdx]))*maxnz[submatIdx] , hackOffsets[submatIdx + 1]);
            for (int idxS = first + tid; idxS < last; idxS += blockDim.x){
                vals[idxS-first] = AS[idxS];
                cols[idxS-first] = JA[idxS];
            }
            __syncthreads();

            for(int t = tid; t<(last-first)/maxnz[submatIdx]*n; t +=blockDim.x){
                for (int j=0; j<maxnz[submatIdx]; j++){
                    res += vals[(t/n)*maxnz[submatIdx] + j] * coeff[cols[(t/n)*maxnz[submatIdx] + j]*n + (t%n)];
                }
                myRes[(hackSize*submatIdx + counter + t/n)*n + (t%n)] = res;
                res=0;
            }
            
            counter += (last-first)/maxnz[submatIdx];
            __syncthreads();
            
        }    
    }
    
}


performance optimized_cuda_h_ellpack_product_in_bis(h_ellpack_matrix_bis host_mat, matrix multivector, matrix* result){

    // Resetto il device
    checkCudaErrors(cudaDeviceReset());

    double* d_as,*d_multivector,* d_result;
    int* d_ja, *d_maxnz, *d_hackOffset;
    
    // Alloco memoria in GPU
    checkCudaErrors(cudaMalloc((void**)&d_as, host_mat.matDim * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_ja, host_mat.matDim * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_maxnz, (host_mat.numMatrix) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_hackOffset, (host_mat.numMatrix+1) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_multivector, multivector.n * multivector.m * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_result, host_mat.m * multivector.n * sizeof(double)));

    // Copia dati sulla GPU
    checkCudaErrors(cudaMemcpy(d_as, host_mat.AS, host_mat.matDim * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_ja, host_mat.JA, host_mat.matDim * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_maxnz, host_mat.maxnz, (host_mat.numMatrix) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_hackOffset, host_mat.hackOffsets, (host_mat.numMatrix+1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_multivector, multivector.coeff, multivector.n * multivector.m * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset(d_result, 0, host_mat.m * multivector.n * sizeof(double)));

    // Assegno numero di blocchi in base al numero di matrici da gestire
    int gridX = host_mat.numMatrix;
    if(gridX > MAX_GRID_SIZE) gridX = MAX_GRID_SIZE;

    dim3 gridSize(gridX);
    dim3 blockSize(1024);

    // Lancio kernel CUDA misurando il tempo impiegato per l'esecuzione 
    float time;
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));
    optimized_cuda_h_ellpack_product<<<gridSize,blockSize>>>(host_mat.m, multivector.n, d_maxnz, d_as, d_ja, d_hackOffset, host_mat.hackSize, host_mat.numMatrix, host_mat.matDim, d_multivector, d_result);
    checkCudaErrors(cudaEventRecord(stop, 0));
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&time, start, stop));
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    
    // Salvo risultato in CPU
    double * resultData =(double*) malloc(host_mat.m * multivector.n * sizeof(double));
    checkCudaErrors(cudaMemcpy(resultData, d_result, host_mat.m * multivector.n * sizeof(double), cudaMemcpyDeviceToHost));
    result->coeff = resultData;

    // Libero memoria in GPU
    checkCudaErrors(cudaFree(d_ja));
    checkCudaErrors(cudaFree(d_as));
    checkCudaErrors(cudaFree(d_multivector));
    checkCudaErrors(cudaFree(d_maxnz));
    checkCudaErrors(cudaFree(d_hackOffset));
    checkCudaErrors(cudaFree(d_result));

    // Trovo Bandwidth effettiva del kernel
    double Br = 8*(host_mat.matDim + multivector.m*multivector.n) + 4*(host_mat.matDim+2*host_mat.numMatrix+1);
    double Bw = (result->m * result->n) * 8;
    double effective_bandwidth = ((Br+Bw)/(pow (10 ,9)))/(time/1000);

    // Restituisco struttura contenente tempo impiegato e bandwidth
    performance perf;
    perf.time = (double)time/1000;
    perf.bandwidth = effective_bandwidth;

    return perf;
}


double optimized_cuda_ellpack_product_in(ellpack_matrix host_mat, matrix vector, matrix* result){

    // Resetto il device
    checkCudaErrors(cudaDeviceReset());

    double* d_as,* d_multivector, * d_result;
    long* d_ja;

    // Alloco memoria in GPU
    checkCudaErrors(cudaMalloc((void**)&d_as, host_mat.m * host.mat.n * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_ja, host_mat.m * host.mat.n * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_multivector, vector.m * vector.n * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_result, host_mat.m * multivector.n * sizeof(double)));

    // Copia dati sulla GPU
    checkCudaErrors(cudaMemcpy(d_as, host_mat.AS, host_mat.m * host.mat.n * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_ja, host_mat.JA, host_mat.m * host.mat.n * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_multivector, multivector.coeff, vector.m * vector.n * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset(d_result, 0, host_mat.m * multivector.n * sizeof(double)));
    
    // Assegno numero di blocchi in modo che ogni blocco gestisca un numero di righe pari a blockX
    int blockX = int(((1024/vector.m)/32)*32), blockY = vector.m;
    if (blockX==0) blockX++;
    int gridX = int(host_mat.m/blockX)+1;

    dim3 gridSize(gridX);
    dim3 blockSize(blockX*blockY);    

    // Lancio kernel CUDA misurando il tempo impiegato per l'esecuzione 
    float time;
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    optimized_cuda_ellpack_product<<<gridSize,blockSize>>>(result->m, result->n, host_mat.maxnz, AS, JA, coeff, cuda_result);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);

    // Salvo risultato in CPU
    double* resultData =(double*) malloc(result->n * result->m * sizeof(double));
    checkCudaErrors(cudaMemcpy(resultData, d_result, result->n * result->m * sizeof(double), cudaMemcpyDeviceToHost));
    result->coeff = resultData;
    
    // Libero memoria in GPU
    checkCudaErrors(cudaFree((void*)JA));
    checkCudaErrors(cudaFree((void*)AS));
    checkCudaErrors(cudaFree((void*)coeff));
    checkCudaErrors(cudaFree((void*)cuda_result));

    return (double)time;
}

#endif