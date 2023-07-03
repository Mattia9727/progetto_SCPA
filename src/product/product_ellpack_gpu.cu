#ifndef _PRODUCTELLPACKGPUH_
#define _PRODUCTELLPACKGPUH_

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "../src/matrices/format/headers/ellpack.h"
#include "../src/matrices/headers/matrix_generator.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>

#define MAX_GRID_SIZE 65536
#define WARP_SIZE 32
#define min(a,b)             \
({                           \
    __typeof__ (a) _a = (a); \
    __typeof__ (b) _b = (b); \
    _a < _b ? _a : _b;       \
})

double* move_AS_on_gpu(ellpack_matrix* mat){
    double* AS;
    checkCudaErrors(cudaMalloc((void**)&AS, sizeof(double) * mat->maxnz));
    checkCudaErrors(cudaMemcpy(AS, mat->AS,  sizeof(double) * mat->maxnz, cudaMemcpyHostToDevice));
    return AS;
}

long* move_JA_on_gpu(ellpack_matrix* mat){
    long* JA;
    checkCudaErrors(cudaMalloc((void**)&JA, sizeof(long) *mat->maxnz));
    checkCudaErrors(cudaMemcpy(JA, mat->JA,  sizeof(long) * mat->maxnz, cudaMemcpyHostToDevice));
    return JA;
}

double* h_move_AS_on_gpu(h_ellpack_matrix* mat){
    double* AS;
    checkCudaErrors(cudaMalloc((void**)&AS, sizeof(double) * mat->matDim));
    checkCudaErrors(cudaMemcpy(AS, mat->AS,  sizeof(double) * mat->matDim, cudaMemcpyHostToDevice));
    return AS;
}

long* h_move_JA_on_gpu(h_ellpack_matrix* mat){
    long* JA;
    checkCudaErrors(cudaMalloc((void**)&JA, sizeof(long) *mat->matDim));
    checkCudaErrors(cudaMemcpy(JA, mat->JA,  sizeof(long) * mat->matDim, cudaMemcpyHostToDevice));
    return JA;
}
long* h_move_maxnz_on_gpu(h_ellpack_matrix* mat){
    long* maxnz;
    checkCudaErrors(cudaMalloc((void**)&maxnz, sizeof(long) *mat->numMatrix));
    checkCudaErrors(cudaMemcpy(maxnz, mat->maxnz,  sizeof(long) * mat->numMatrix, cudaMemcpyHostToDevice));
    return maxnz;
}

long* h_move_hackOffsets_on_gpu(h_ellpack_matrix* mat){
    long* hackOffsets;
    checkCudaErrors(cudaMalloc((void**)&hackOffsets, sizeof(long) *(mat->numMatrix+1)));
    checkCudaErrors(cudaMemcpy(hackOffsets, mat->hackOffsets,  sizeof(long) * (mat->numMatrix+1), cudaMemcpyHostToDevice));
    return hackOffsets;
}

double* move_matrix_coeff_on_gpu(matrix* vector){
    double* coeff;
    checkCudaErrors(cudaMalloc((void**)&coeff, sizeof(double) * vector->n * vector->m));
    checkCudaErrors(cudaMemcpy(coeff, vector->coeff,  sizeof(double) * vector->n * vector->m, cudaMemcpyHostToDevice));
    return coeff;
}

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


__global__ void optimized_cuda_h_ellpack_product(int m, int n, long* maxnz, double* AS, long* JA, long* hackOffsets, long hackSize, long numMatrix, long matDim, double* coeff, double* myRes){
    printf("ciao belli\n");
    long idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx <= m*n){
        long AS_JA_row = idx/n;                     //Indice di riga della matrice risultante (con cui mi scorro AS e JA)
        long multiv_col = idx%n;                    //Indice di colonna della matrice risultante (con cui mi scorro il multivettore)
        double res=0;
        long i;
        long rowCount=0;
        for (i=0; i<matDim; i++){
            rowCount+=hackSize*n;
            if (idx < rowCount) break;              //Voglio recuperare l'indice della sottomatrice giusta in modo da avere il relativo maxnz e hackOffset.
        }
        long subrow = AS_JA_row - hackSize*i;
        
        for (long j=0; j<maxnz[i]; j++){
            //if (AS[hackOffsets[i] + subrow*maxnz[i] + j] == 0) break;
            res+=AS[hackOffsets[i] + subrow*maxnz[i] + j] * coeff[multiv_col*m + JA[hackOffsets[i] + subrow*maxnz[i] + j]];
        }
        myRes[AS_JA_row * n + multiv_col] = res;

    }
}
__device__ double warp_reduce_2(double val){
    for(int offset=warpSize/2; offset>0; offset/=2){
        val += __shfl_down_sync(0xffffffff,val, offset);
    }
    return val;
}

__global__ void optimized_cuda_h_ellpack_product_bis(int m, int nCols, long* maxnz, double* AS, long* JA, long* hackOffsets, int hackSize, long numMatrix, double* multivector, double* myRes){
   long idxBlock = blockIdx.x; // Indice del blocco che determina la sottomatrice da fare
   long warpId = threadIdx.y; //Indice del warp
   long tid = threadIdx.x; //Indice del thread nel warp
   long start,end;
   double val;
   long col;
   double sum[2048] = {0};
   for(int subMat = idxBlock; subMat < numMatrix; subMat += MAX_GRID_SIZE){
    // Devo identificare la sottomatrice
    start = hackOffsets[subMat];
    end = hackOffsets[subMat+1];
    //if(m - idxBlock*hackSize - warpId > 0){
        for(int idx = start + warpId*maxnz[idxBlock]; idx < start + (warpId+1)*maxnz[idxBlock]; idx +=WARP_SIZE){
            val = AS[hackOffsets[idxBlock]+warpId*maxnz[idxBlock] + idx];
            col = JA[hackOffsets[idxBlock]+idxBlock*warpId*maxnz[idxBlock] + idx];
            for(int col_m = 0; col_m < nCols; col_m++){
                sum[warpId*nCols + col_m] += val * multivector[col*nCols + col_m];
            }
        }
    //}
    for(int j = 0; j < nCols; j++){
        sum[warpId*nCols +j] = warp_reduce_2(sum[warpId*nCols +j]);
        if(tid == 0){
            myRes[(idxBlock*blockDim.x + warpId)*nCols+j] = sum[warpId*nCols +j];
            sum[warpId*nCols +j] = 0;
        }
    }
   }

}


double optimized_cuda_ellpack_product_in(ellpack_matrix host_mat, matrix vector, matrix* result){

    printf("maxnz=%d\n", host_mat.maxnz);

    double* AS = move_AS_on_gpu(&host_mat);
    long* JA = move_JA_on_gpu(&host_mat);
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
    checkCudaErrors(cudaGetLastError());
    cudaEventRecord(stop, 0);// 0 - the default stream
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("%f\n",time);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    result->coeff = myRes;
    checkCudaErrors(cudaFree((void*)JA));
    checkCudaErrors(cudaFree((void*)AS));
    checkCudaErrors(cudaFree((void*)coeff));

    cudaMemcpy(result->coeff, cuda_result, sizeof(double) * result->n * result->m, cudaMemcpyDeviceToHost);
    
    checkCudaErrors(cudaFree((void*)cuda_result));
    return (double)time;
}

double optimized_cuda_h_ellpack_product_in(h_ellpack_matrix host_mat, matrix vector, matrix* result){
    //printf("maxnz: ");
    //for(int i=0; i<host_mat.numMatrix;i++) printf("%ld ", host_mat.maxnz[i]);
    //printf("\n");

    // resetto il device
    checkCudaErrors(cudaDeviceReset());

    printf("m=%d, n=%d\n",result->m,result->n);

    double* AS = h_move_AS_on_gpu(&host_mat);
    long* JA = h_move_JA_on_gpu(&host_mat);
    long* maxnz = h_move_maxnz_on_gpu(&host_mat);
    long* hackOffsets = h_move_hackOffsets_on_gpu(&host_mat);
    double* coeff = move_matrix_coeff_on_gpu(&vector);

    double* cuda_result;
    checkCudaErrors(cudaMalloc((void**)&cuda_result, sizeof(double) * result->m * result->n));
    
    int blockX = int(((1024/vector.m)/32)*32), blockY = vector.m;  //gestione warp. Ogni blocco gestisce un numero di righe pari a blockX
    if (blockX==0) blockX++;
    int gridX = int(host_mat.m/blockX)+1;

    //printf("blockx=%d, blocky=%d, gridx=%d\n", blockX, blockY, gridX);

    //dim3 gridSize(gridX);
    //dim3 blockSize(blockX*blockY);    
    printf("%d\n",host_mat.numMatrix);
    dim3 gridSize((int)host_mat.numMatrix);               //NUMERO DI BLOCCHI IN UNA GRID
    dim3 blockSize(32,32);          //NUMERO DI THREAD IN UN BLOCCO

    float time;
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));// 0 - the default stream
    optimized_cuda_h_ellpack_product_bis<<<gridSize,blockSize>>>(result->m, result->n, maxnz, AS, JA, hackOffsets, (int)host_mat.hackSize, host_mat.numMatrix, coeff, cuda_result);
    checkCudaErrors(cudaGetLastError());
    checkCudaErrors(cudaDeviceSynchronize());
    checkCudaErrors(cudaEventRecord(stop, 0));// 0 - the default stream
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&time, start, stop));
    printf("time: %f\n\n",time);
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    checkCudaErrors(cudaFree(JA));
    checkCudaErrors(cudaFree(AS));
    checkCudaErrors(cudaFree(coeff));
    checkCudaErrors(cudaFree(maxnz));
    checkCudaErrors(cudaFree(hackOffsets));

    cudaMemcpy(result->coeff, cuda_result, sizeof(double) * result->n * result->m, cudaMemcpyDeviceToHost);
    
    checkCudaErrors(cudaFree((void*)cuda_result));
    return (double)time/1000;
}
#endif