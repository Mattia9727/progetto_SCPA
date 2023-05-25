#include "headers/product_csr.h"
#include <cuda_runtime.h>

#define SHARED_MEMORY_DIM 49152

void calculate_rows_block(int totalRows, int* irp, int* rowBlocks){
    rowBlocks[0] = 0; 
    int sum = 0, last_i= 0, ctr=1; 
    for(int i = 1; i < totalRows; i++){
        sum += irp[i]-irp[i-1];
        if(sum == SHARED_MEMORY_DIM/sizeof(double)){
            last_i = i;
            rowBlocks[ctr++] = i;
            sum = 0;

        }
        else if( sum > SHARED_MEMORY_DIM/sizeof(double)){
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
    rowBlocks[ctr++] = totalRows;
}                                           

__global__ void hello(char *a, int *b) 
{
	a[threadIdx.x] += b[threadIdx.x];
}
 



double calcola_prodotto_csr_cuda(csr_matrix mat, matrix multivector, matrix* result){
    int N = 16;
    int blocksize = 16;
    char a[N] = "Hello \0\0\0\0\0\0";
	int b[N] = {15, 10, 6, 0, -11, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
 
	char *ad;
	int *bd;
	const int csize = N*sizeof(char);
	const int isize = N*sizeof(int);
 
	printf("%s", a);
 
	cudaMalloc( (void**)&ad, csize ); 
	cudaMalloc( (void**)&bd, isize ); 
	cudaMemcpy( ad, a, csize, cudaMemcpyHostToDevice ); 
	cudaMemcpy( bd, b, isize, cudaMemcpyHostToDevice ); 
	
	dim3 dimBlock( blocksize, 1 );
	dim3 dimGrid( 1, 1 );
	hello<<<dimGrid, dimBlock>>>(ad, bd);
	cudaMemcpy( a, ad, csize, cudaMemcpyDeviceToHost ); 
	cudaFree( ad );
	cudaFree( bd );
	
	printf("%s\n", a);
	return EXIT_SUCCESS;
}