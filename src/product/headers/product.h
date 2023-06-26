#ifndef _PRODUCTH_
#define _PRODUCTH_
#include "../../matrices/format/headers/csr.h"


void prepara_risultato(int m, int n, matrix* result);
void prepara_risultato_cuda(int m, int n, matrix* result);

void free_matrix(matrix* result);
void free_matrix_cuda(matrix* result);

double calcola_prodotto_seriale(csr_matrix csrMatrix, matrix vector, matrix* result);
double check_result(matrix m1, matrix m2);

#endif