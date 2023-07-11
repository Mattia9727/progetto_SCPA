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


__global__ void optimized_cuda_h_ellpack_product(int m, int n, int* maxnz, double* AS, int* JA, int* hackOffsets, int hackSize, int numMatrix, int matDim, double* coeff, double* myRes){
    int bid = blockIdx.x;
    int tid = threadIdx.x;
    
    __shared__ double vals[4096];
    __shared__ int cols[4096];

    int first = 0;
    int last = 0;
    double res=0;
    for (int submatIdx = bid; submatIdx<numMatrix; submatIdx += gridDim.x){
        if(maxnz[submatIdx] < 4096){
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
        }else{
            for(int t = tid; t<hackSize*n; t +=blockDim.x){
                for (int j=0; j<maxnz[submatIdx]; j++){
                    res += AS[hackOffsets[submatIdx]+(t/n)*maxnz[submatIdx] + j] * coeff[JA[hackOffsets[submatIdx]+(t/n)*maxnz[submatIdx] + j]*n + (t%n)];
                }
                myRes[(hackSize*submatIdx + t/n)*n + (t%n)] = res;
                res=0;
            }
            __syncthreads();
        }
        
    }
    
}

__device__ double warp_reduce_2(double val){
    for(int offset=warpSize/2; offset>0; offset/=2){
        val += __shfl_down_sync(0xffffffff,val, offset);
    }
    return val;
}
/*
__global__ void optimized_cuda_h_ellpack_product_2(int nRows, int nCols, int* maxnz, double* AS, int* JA, int* hackOffsets, int hackSize, int numMatrix, double* multivector, double* myRes){
    int idxBlock = blockIdx.x; // Indice del blocco che determina la sottomatrice da fare
    int warpId = threadIdx.y; //Indice del warp
    int tid = threadIdx.x; //Indice del thread nel warp
    int start,end;
    double val;
    int col,maxnzBlock,col_mul,col_incr = 1;
    int nReplic = 0, state = 0;
    __shared__ int sharedState;
    __shared__ double vals[3392];
    __shared__ int cols[3392]; 
    __shared__ double sum[1024];
    // Ciclo nel caso in cui il num di blocchi fosse minore del num di righe
    for(int subMat = idxBlock; subMat < numMatrix; subMat += gridDim.x){
        // Devo identificare la sottomatrice    
        start = hackOffsets[subMat];
        end = hackOffsets[subMat+1];
        maxnzBlock = maxnz[subMat];

        if(maxnzBlock <= 106){
            sharedState = 0;
            if(maxnzBlock <= 32 && maxnzBlock >0){
                nReplic = (int)32/maxnzBlock;
                col_incr = nReplic;
            }
        }else sharedState = -1;

        for(int col_m = 0; col_m < nCols; col_m+=col_incr){
            //questo  ciclo mi scorre nelle righe della sottomatrice
            for(int warp = warpId; warp < hackSize; warp += WARP_SIZE){
                if(nRows - subMat*hackSize - warp > 0){
                    sum[warpId*32 + tid] = 0.0;
                    //questo ciclo mi scorre negli elementi della riga
                    if(sharedState == 0 || nReplic <= 0){
                        for(int idx = start + warp*maxnzBlock+tid; idx < start + (warp+1)*maxnzBlock; idx +=WARP_SIZE){
                            if(sharedState == -1){
                                val = AS[idx];
                                col = JA[idx];
                                sum[warp*32 + tid] += val * multivector[col*nCols + col_m];
                            }                                
                            else if(sharedState == 0){
                                val = AS[idx];
                                col = JA[idx];
                                vals[warp*106 +idx-start-warp*maxnzBlock] = val;
                                cols[warp*106 +idx-start-warp*maxnzBlock] = col;
                                if(nReplic <= 0)
                                    sum[warp*32 + tid] += val * multivector[col*nCols + col_m];
                                state = 1;
                            }else{
                                sum[warp*32 + tid] += vals[warp*106+idx-start-warp*maxnzBlock] * multivector[cols[warp*106+idx-start-warp*maxnzBlock]*nCols + col_m];
                            }
                        }
                    }
                    if(sharedState == 0){
                        __syncthreads();
                        if(state == 1) sharedState = 1;
                        __syncthreads();
                    }
                    if(sharedState == 1 && nReplic > 0){
                        //shared state pronta 
                        if(tid < nReplic*maxnzBlock){
                            col_mul = tid/maxnzBlock + col_m;
                            if(col_mul < nCols){
                                sum[warp*32 + tid] += vals[warp*106+tid%maxnzBlock] * multivector[cols[warp*106+tid%maxnzBlock]*nCols + col_mul];
                            }
                                
                        }
                    }  
                }
            }
            
            __syncthreads();
            
            if(tid == 0 && nReplic <= 0){
                for(int i = 1; i < 32; i++){
                    sum[warpId*32] += sum[warpId*32+i];
                }
                myRes[(subMat*hackSize + warpId)*nCols+col_m] = sum[warpId*32];   
            } else if(nReplic > 0 && tid < nReplic){
                for(int i = 1; i < maxnzBlock; i++){
                    sum[warpId*32+tid*maxnzBlock] += sum[warpId*32 + tid*maxnzBlock + i];
                }
                if(col_m + tid < nCols){
                    myRes[(subMat*hackSize + warpId)*nCols+col_m + tid] = sum[warpId*32+tid*maxnzBlock];
                }
                
            }
            __syncthreads();
            
        }
        
    }

}
*/
performance optimized_cuda_h_ellpack_product_in_bis(h_ellpack_matrix_bis host_mat, matrix multivector, matrix* result){

    // resetto il device
    checkCudaErrors(cudaDeviceReset());

    double* d_as,*d_multivector,* d_result;
    int* d_ja, *d_maxnz, *d_hackOffset;
    
    // Alloco
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

    int gridX = host_mat.numMatrix;
    if(gridX > MAX_GRID_SIZE) gridX = MAX_GRID_SIZE;

    dim3 gridSize(gridX);               //NUMERO DI BLOCCHI IN UNA GRID
    dim3 blockSize(1024);              //NUMERO DI THREAD IN UN BLOCCO

    float time;
    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));// 0 - the default stream
    optimized_cuda_h_ellpack_product<<<gridSize,blockSize>>>(host_mat.m, multivector.n, d_maxnz, d_as, d_ja, d_hackOffset, host_mat.hackSize, host_mat.numMatrix, host_mat.matDim, d_multivector, d_result);
        //optimized_cuda_h_ellpack_product<<<gridSize,blockSize>>>(result->m, result->n, maxnz, AS, JA, hackOffsets, host_mat.hackSize, host_mat.numMatrix, host_mat.matDim, coeff, cuda_result);
    checkCudaErrors(cudaEventRecord(stop, 0));// 0 - the default stream
    checkCudaErrors(cudaEventSynchronize(stop));
    checkCudaErrors(cudaEventElapsedTime(&time, start, stop));
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    checkCudaErrors(cudaFree(d_ja));
    checkCudaErrors(cudaFree(d_as));
    checkCudaErrors(cudaFree(d_multivector));
    checkCudaErrors(cudaFree(d_maxnz));
    checkCudaErrors(cudaFree(d_hackOffset));

    double * resultData =(double*) malloc(host_mat.m * multivector.n * sizeof(double));
    
    checkCudaErrors(cudaMemcpy(resultData, d_result, host_mat.m * multivector.n * sizeof(double), cudaMemcpyDeviceToHost));
    
    result->coeff = resultData;
    
    checkCudaErrors(cudaFree(d_result));
    performance perf;
    perf.time = (double)time/1000;
    double Br = 8*(host_mat.matDim/perf.time  + (multivector.m*multivector.n)/perf.time) + 4*(6/perf.time + host_mat.matDim/perf.time+(2*host_mat.numMatrix/perf.time));
    double Bw = ((result->m * result->n)/perf.time) * 8;
    double effective_bandwidth = ((Br+Bw)/(pow (10 ,9)));
    perf.bandwidth = effective_bandwidth;
    return perf;
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
    //optimized_cuda_h_ellpack_product<<<gridSize,blockSize>>>(result->m, result->n, maxnz, AS, JA, hackOffsets, host_mat.hackSize, host_mat.numMatrix, host_mat.matDim, coeff, cuda_result);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);// 0 - the default stream
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("time: %f\n\n",time);
    cudaEventDestroy(start);
    cudaEventDestroy(stop);
    result->coeff = myRes;
    checkCudaErrors(cudaFree(JA));
    checkCudaErrors(cudaFree(AS));
    checkCudaErrors(cudaFree(coeff));
    checkCudaErrors(cudaFree(maxnz));
    checkCudaErrors(cudaFree(hackOffsets));

    cudaMemcpy(result->coeff, cuda_result, sizeof(double) * result->n * result->m, cudaMemcpyDeviceToHost);

    checkCudaErrors(cudaFree((void*)cuda_result));
    double Br = (host_mat.matDim + vector.m*vector.n) * 8;
    double Bw = (result->m * result->n) * 8;
    double effective_bandwidth = ((Br+Bw)/(pow (10 ,9)))/(time/1000);
    printf("EB: %f\n",effective_bandwidth);
    return (double)time/1000;
}
#endif