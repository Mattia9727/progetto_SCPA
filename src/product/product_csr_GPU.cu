#include "../matrices/format/headers/csr.h"

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

void calcola_prodotto_csr_cuda(csr_matrix mat, matrix multivector, matrix* result){
    
}