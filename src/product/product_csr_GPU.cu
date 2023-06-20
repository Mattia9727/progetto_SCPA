#include "headers/product_csr.h"
#include <cuda_runtime.h>
#include <helper_cuda.h>
#include <iostream>

#define MAX_NNZ_PER_WG 6144
#define MAX_BLOCK_THREADS 1024
#define WARP_SIZE 32

__device__ double warp_reduce(double val){
    for(int offset=warpSize/2; offset>0; offset/=2){
        val += __shfl_down_sync(0xffffffff,val, offset);
    }
    return val;
}
 
int* calculate_rows_block(int totalRows, int* irp, int* numBlocks){
    int* rowBlocks = (int*)malloc(21000*sizeof(int));
    rowBlocks[0] = 0; 
    int sum = 0, last_i= 0, ctr=1; 
    for(int i = 1; i < totalRows; i++){
        sum += irp[i]-irp[i-1];
        if(sum == MAX_NNZ_PER_WG){
            last_i = i;
            rowBlocks[ctr++] = i;
            sum = 0;

        }
        else if( sum > MAX_NNZ_PER_WG){
            if(i - last_i > 1){
                rowBlocks[ctr++] = i-1;
                i--;
            }
            else if(i - last_i == 1){
                rowBlocks[ctr++] = i;
            }
            last_i = i; 
            sum = 0;
        } 
    }
    
    //printf("%d %d\n",ctr,totalRows);
    *numBlocks = ctr;
    rowBlocks[ctr++] = totalRows;
    return rowBlocks;
} 

__global__ void gpuMultiplication(double* as, int* ja, int* irp, double* multivector, int m, int n, int col_multivector, double* resultData){
    //CSR-Vector
    extern __shared__ double vals[];
 
    int threadId = blockIdx.x * blockDim.x + threadIdx.x; // global thread index
    int warpId = threadId / 32; // Global warp index
    int lane = threadId &(32-1); // thread index within the warp
    //one warp per row
    int row = warpId;
    int startRow = irp[row];
    int stopRow = irp[row+1];
    int col;
    double val;
   
    if(row < m){
        for(int i = startRow + lane; i < stopRow; i +=32){
            val = as[i];
            col = ja[i];
            
            for(int col_m = 0; col_m < col_multivector; col_m++){
                vals[col_m] = 0;
                vals[col_m] += val * multivector[col*col_multivector + col_m];
                vals[col_m] = warp_reduce(vals[col_m]);
                //first thread writes the result 
                if(lane == 0)
                    resultData[row* col_multivector + col_m] += vals[col_m];
            }    
        }
    }
      
}

__global__ void sparseDenseMatrixMul(double* as, int* ja, int* irp, int m, double* multivector, int n, int q, double* resultData)
{
    int row_mat = blockIdx.x * blockDim.x + threadIdx.x;
    int col_vec = threadIdx.y;
    

    if (row_mat < m && col_vec < q)
    {

        int start = irp[row_mat];
        int end = irp[row_mat+1];
        double sum = 0;
        
        
        for (int i = start; i < end; i++)
        {
            sum += as[i] * multivector[ja[i] * q + col_vec];
        }
        
        resultData[row_mat * q + col_vec] = sum;
    }
}


__device__ int prev_power_of_2(int x){
    x = x | (x >> 1);
    x = x | (x >> 2);
    x = x | (x >> 4);
    x = x | (x >> 8);
    x = x | (x >> 16);
    return x - (x >> 1);
}

__global__ void csrAdaptiveMult(double* as, int* ja, int* irp, double* multivector, int m, int n, int col_multivector, int* rowBlocks, double* resultData){
    __shared__ double vals[MAX_NNZ_PER_WG];
    
    int startRow = rowBlocks[blockIdx.x];
    int stopRow = rowBlocks[blockIdx.x+1];
    long int numRows = stopRow - startRow;
    int nnz = irp[stopRow]-irp[startRow];
    if (numRows > 1){
        //CSR-Stream
        //printf("csr stream\n");
        int tid = threadIdx.x; // indice del thread nel blocco
        int localCol;
        for(int j = 0; j < col_multivector; j++){
            for(int i = tid; i < nnz; i+= blockDim.x){ 
                localCol = irp[startRow]+i;
                vals[i] = as[localCol];
                vals[i] *= multivector[ja[localCol]*col_multivector+j];
            }
            int firstCol = irp[startRow];
            int localRow = startRow + tid;
            
            __syncthreads();
            while(localRow < stopRow){
                double temp = 0; 
                for(int i = irp[localRow]-firstCol; i < irp[localRow+1]-firstCol; i++){
                    temp += vals[i];
                }
                resultData[localRow*col_multivector +j] = temp;
                localRow += blockDim.x;
            }
            __syncthreads();
        }
    }else {
        //CSR-Vector
        int threadId = threadIdx.x; // global thread index
        //printf("csr vector\n");
        int warpId = threadId / 32; // Global warp index
        int lane = threadId &(32-1); // thread index within the warp
        //one warp per row
        double val; 
        int col;
        double sum[64] = {0};
        //Questo blocco fa solo questa riga
        if(warpId == 0){            
            for(int i = irp[startRow] + lane; i < irp[startRow+1]; i +=32){
                val = as[i];
                col = ja[i];
                for(int j = 0; j < col_multivector; j++){
                    sum[j] += val*multivector[col*col_multivector + j];
                }    
            }
        }
        for(int i = 0; i < col_multivector; i++){
            sum[i] = warp_reduce(sum[i]);
            if(lane == 0 && warpId == 0){
                
                resultData[startRow*col_multivector + i] = sum[i];   
            }    
        }       
    }
}

double calcola_prodotto_csr_cuda(csr_matrix mat, matrix multivector, matrix* result){
    
    int num_blocks;
    int* res = calculate_rows_block(mat.m, mat.irp,&num_blocks);

    double* d_as;
    int* d_ja;
    int* d_irp, *d_rowBlocks;
    double* d_multivector, *d_multivector2;
    double* d_result, *d_result2;

    // resetto il device
    checkCudaErrors(cudaDeviceReset());

    // Definizione della griglia, del blocco e della shared memory

    int blockX = MAX_BLOCK_THREADS;
    int gridX = num_blocks;
    /*
    if(mat.m * WARP_SIZE > MAX_BLOCK_THREADS){
        gridX = (mat.m / (MAX_BLOCK_THREADS/WARP_SIZE)) +1;
    }
    if(gridX < num_blocks-1) gridX = num_blocks-1;
    */
    dim3 blockSize(blockX);
    dim3 gridSize(gridX);
    //sharedMemorySize = sizeof(double)*MAX_NNZ_PER_WG;
    // Alloco
    checkCudaErrors(cudaMalloc((void**)&d_as, mat.nz * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_ja, mat.nz * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_irp, (mat.m + 1) * sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_rowBlocks, (num_blocks+1) * sizeof(int)));

    checkCudaErrors(cudaMalloc((void**)&d_multivector, multivector.m * multivector.n * sizeof(double)));
    checkCudaErrors(cudaMalloc((void**)&d_result, mat.m * multivector.n * sizeof(double)));


    // Copia dati sulla GPU
    checkCudaErrors(cudaMemcpy(d_as, mat.as, mat.nz * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_ja, mat.ja, mat.nz * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_irp, mat.irp, (mat.m + 1) * sizeof(int), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemcpy(d_rowBlocks, res, (num_blocks+1) * sizeof(int), cudaMemcpyHostToDevice));
    
    checkCudaErrors(cudaMemcpy(d_multivector, multivector.coeff, multivector.m * multivector.n * sizeof(double), cudaMemcpyHostToDevice));
    checkCudaErrors(cudaMemset(d_result, 0, mat.m * multivector.n * sizeof(double)));


    cudaEvent_t start, stop;
    checkCudaErrors(cudaEventCreate(&start));
    checkCudaErrors(cudaEventCreate(&stop));
    checkCudaErrors(cudaEventRecord(start, 0));
    
    // Esecuzione del kernel
    csrAdaptiveMult<<<gridSize, blockSize>>>(d_as, d_ja, d_irp, d_multivector, mat.m,mat.n, multivector.n, d_rowBlocks,d_result);

    checkCudaErrors(cudaEventRecord(stop, 0));
    checkCudaErrors(cudaEventSynchronize(stop));
    //sparseDenseMatrixMul<<<gridSize, blockSize>>>(d_as, d_ja, d_irp, mat.m, d_multivector, mat.n, multivector.n, d_result);
    float time;
    checkCudaErrors(cudaEventElapsedTime(&time, start, stop));
    checkCudaErrors(cudaEventDestroy(start));
    checkCudaErrors(cudaEventDestroy(stop));
    // Copia risultato dalla GPU alla CPU
    double * resultData =(double*) malloc(mat.m * multivector.n * sizeof(double));
    
    checkCudaErrors(cudaMemcpy(resultData, d_result, mat.m * multivector.n * sizeof(double), cudaMemcpyDeviceToHost));
    
    result->coeff = resultData;
    // Deallocazione memoria
    checkCudaErrors(cudaFree(d_ja));
    checkCudaErrors(cudaFree(d_as));
    checkCudaErrors(cudaFree(d_irp));
    checkCudaErrors(cudaFree(d_multivector));
    checkCudaErrors(cudaFree(d_result));
    checkCudaErrors(cudaFree(d_rowBlocks));

    return (double)time/1000;
    
}