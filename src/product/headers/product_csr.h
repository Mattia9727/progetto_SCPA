
#ifndef _PRODUCTCSROPENMPH_
#define _PRODUCTCSROPENMPH_ 

#include "product.h"

performance calcola_prodotto_per_righe_csr_openmp(csr_matrix csrMatrix, matrix multivector,matrix* result, int nThreads);

double calcola_prodotto_per_righe_csr_openmp_bis(csr_matrix csrMatrix, matrix multivector,matrix* result, int nThreads);
double calcola_prodotto_per_righe_csr_openmp_trasposto(csr_matrix csrMatrix, matrix multivector_trasposto,matrix* result, int nThreads);
performance calcola_prodotto_csr_cuda(csr_matrix mat, matrix multivector, matrix* result);    
#endif