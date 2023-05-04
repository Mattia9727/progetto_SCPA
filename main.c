#include <malloc.h>
#include <stdlib.h>
#include "product/product_csr_openmp.h"
#include <time.h>

#define NFILES 8

void main(){
    coo_matrix mat;
    csr_matrix converted_matrix;
    matrix multivector, result;
    clock_t begin, end;
    int col_multivector[7] = {3,4,8,12,16,32,64}; 

    char* matFiles[] = {
        "test_matrices/matrix_files/cage4.mtx",
        "test_matrices/matrix_files/adder_dcop_32.mtx",
        "test_matrices/matrix_files/af23560.mtx",
        "test_matrices/matrix_files/bcsstk17.mtx",
        "test_matrices/matrix_files/cant.mtx",
        "test_matrices/matrix_files/cavity10.mtx",
        "test_matrices/matrix_files/cop20k_A.mtx",
        "test_matrices/matrix_files/Cube_Coup_dt0.mtx",
        "test_matrices/matrix_files/lung2.mtx",
        "test_matrices/matrix_files/mac_econ_fwd500.mtx",
        "test_matrices/matrix_files/mcfe.mtx",
        "test_matrices/matrix_files/mhd4800a.mtx",
        "test_matrices/matrix_files/mhda416.mtx",
        "test_matrices/matrix_files/ML_Laplace.mtx",
        "test_matrices/matrix_files/olafu.mtx",
        "test_matrices/matrix_files/olm1000.mtx",
        "test_matrices/matrix_files/PRO2R.mtx",
        "test_matrices/matrix_files/west2021.mtx",
        "test_matrices/matrix_files/rdist2.mtx",
        "test_matrices/matrix_files/raefsky2.mtx",
    }; 
    
    for(int i = 0; i < NFILES; i++){
        printf("%s\n",matFiles[i]);
        mat = getMatrix(matFiles[i]);
        converted_matrix = convertToCsrFromCoo(mat);
        //stampaMatriceCsr(converted_matrix);

        for(int i = 0; i < 7; i++){
            printf("COLONNE MULTIVETTORE %d\n",col_multivector[i]);
            multivector = GenerateMultivector(mat.n,col_multivector[i]); 
            begin = clock();
            result = prepara_risultato(converted_matrix.m, multivector.n);
            calcola_prodotto_seriale(converted_matrix,multivector, &result);
            end = clock();
            printf("\tserial %f\n",(double)(end - begin) / CLOCKS_PER_SEC);
            // Stampa della matrice
            //stampaMatrice(converted_matrix.m, vector.n,result);
            begin = clock();
            calcola_prodotto_per_righe_csr_openmp(converted_matrix,multivector, &result);
            end = clock();
            printf("\tparall %f\n", (double)(end - begin) / CLOCKS_PER_SEC);
            // Stampa della matrice
            //stampaMatrice(converted_matrix.m, vector.n,result);
    
        }
    }
    /*
    mat = getMatrix("test_matrices/matrix_files/prova_simm.mtx");
    
    matrix new_matrix;
    new_matrix.m = mat.m;
    new_matrix.n = mat.n;
    new_matrix.coeff = (float **)malloc(sizeof(float*)*mat.m);
    for(int i = 0; i < mat.m; i++){
        new_matrix.coeff[i] = (float *) calloc(mat.n, sizeof(float));     // allocazione della memoria per la riga i-esima
        if (new_matrix.coeff[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            exit(0);
        }
    }
    if (new_matrix.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }
    for(int i = 0; i < mat.nz; i++){
        new_matrix.coeff[mat.rows[i]][mat.cols[i]] = (float)(mat.values[i]);
    }
    stampaMatrice(new_matrix);
    */
}