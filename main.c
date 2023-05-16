#include <malloc.h>
#include <stdlib.h>
#include "product/product_csr_openmp.h"
#include <time.h>
#include <math.h>
#include "matrices/format/ellpack.h"
#include "product/ellpack_product.h"

#define NFILES 30
#define NREPETITIONS 50
#define FILEHEADER "MatrixName, NRows, NCol, NZ, k, MFLOPS\n"

void main(){
    coo_matrix mat;
    csr_matrix converted_csr_matrix;
    ellpack_matrix converted_ellpack_matrix;
    matrix multivector, resultCSR, resultELLPACK;
    clock_t begin, end;
    double timeSumSer = 0.0, timeSumCsrOmp = 0.0, timeSumEllpackOmp = 0.0;
    int col_multivector[7] = {3,4,8,12,16,32,64}; 

    FILE *resultsSer = fopen("measurements/results/serial.csv", "w+");
    FILE *resultsCsrOmp = fopen("measurements/results/crsomp.csv", "w+");
    FILE *resultsEllpackOmp = fopen("measurements/results/ellpackomp.csv", "w+");

    if(resultsSer == NULL || resultsCsrOmp == NULL || resultsEllpackOmp == NULL )
        printf("errore apertura file\n");

    char* matFiles[] = {
        "matrices/matrix_files/cage4.mtx",
        "matrices/matrix_files/olm1000.mtx",
        "matrices/matrix_files/west2021.mtx",
        "matrices/matrix_files/mhda416.mtx",
        "matrices/matrix_files/adder_dcop_32.mtx",
        "matrices/matrix_files/mcfe.mtx",
        "matrices/matrix_files/rdist2.mtx",
        "matrices/matrix_files/cavity10.mtx",
        "matrices/matrix_files/mhd4800a.mtx",
        "matrices/matrix_files/bcsstk17.mtx",
        "matrices/matrix_files/raefsky2.mtx",
        "matrices/matrix_files/thermal1.mtx",
        "matrices/matrix_files/af23560.mtx",
        "matrices/matrix_files/thermomech_TK.mtx",
        "matrices/matrix_files/olafu.mtx",
        "matrices/matrix_files/FEM_3D_thermal1.mtx",
        "matrices/matrix_files/lung2.mtx",
        "matrices/matrix_files/dc1.mtx",
        "matrices/matrix_files/amazon0302.mtx",
        "matrices/matrix_files/roadNet-PA.mtx",
        "matrices/matrix_files/cop20k_A.mtx",
        "matrices/matrix_files/mac_econ_fwd500.mtx",
        "matrices/matrix_files/cant.mtx",
        "matrices/matrix_files/webbase-1M.mtx",
        "matrices/matrix_files/thermal2.mtx",
        "matrices/matrix_files/nlpkkt80.mtx",
        "matrices/matrix_files/PR02R.mtx",
        "matrices/matrix_files/af_1_k101.mtx",
        "matrices/matrix_files/ML_Laplace.mtx",
        "matrices/matrix_files/Cube_Coup_dt0.mtx",
    }; 

    fprintf(resultsSer,"%s",FILEHEADER);
    fprintf(resultsCsrOmp,"%s",FILEHEADER);
    fprintf(resultsEllpackOmp,"%s",FILEHEADER);

    for(int i = 0; i < NFILES; i++){
        printf("%s\n",matFiles[i]);
        mat = getMatrix(matFiles[i]);
        converted_csr_matrix = convertToCsrFromCoo(mat);
        converted_ellpack_matrix = ConvertCOOToELLPACK(mat);
        GenerateMultivector(&multivector, mat.n, 64);
        for(int j = 0; j < 7; j++){
            multivector.n=col_multivector[j];
            printf("    Dim: %dx%d, Colonne %d\n",mat.m, mat.n, col_multivector[j]);
            fflush(stdout);
            resultCSR = prepara_risultato(converted_csr_matrix.m, multivector.n);
            resultELLPACK = prepara_risultato(converted_ellpack_matrix.m, multivector.n);
            for(int k = 0; k < NREPETITIONS; k++){
                begin = clock();
                calcola_prodotto_seriale(converted_csr_matrix,multivector, &resultCSR);
                end = clock();
                timeSumSer += (double)(end - begin) / CLOCKS_PER_SEC;
                begin = clock();
                calcola_prodotto_per_righe_csr_openmp(converted_csr_matrix,multivector, &resultCSR);
                end = clock();
                timeSumCsrOmp += (double)(end - begin) / CLOCKS_PER_SEC;
                begin = clock();
                OmpELLPACKProduct(converted_ellpack_matrix, multivector, &resultELLPACK);
                end = clock();
                timeSumEllpackOmp += (double)(end - begin) / CLOCKS_PER_SEC;
            }
            fprintf(resultsSer,"%s, %d, %d, %d, %d, %f\n",matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], ((2*col_multivector[j]/pow(10,6))*converted_csr_matrix.nz)/(timeSumSer/NREPETITIONS));
            fprintf(resultsCsrOmp,"%s, %d, %d, %d, %d, %f\n",matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], ((2*col_multivector[j]/pow(10,6))*converted_csr_matrix.nz)/(timeSumCsrOmp/NREPETITIONS));
            fprintf(resultsEllpackOmp,"%s, %d, %d, %d, %d, %f\n",matFiles[i],converted_ellpack_matrix.m, converted_ellpack_matrix.n, converted_ellpack_matrix.maxnz, col_multivector[j], ((2*col_multivector[j]/pow(10,6))*converted_csr_matrix.nz)/(timeSumEllpackOmp/NREPETITIONS));

            free_matrix(&resultCSR);
            free_matrix(&resultELLPACK);
        }
        free_matrix(&multivector);
        FreeCsrMatrix(&converted_csr_matrix);
        FreeEllpackMatrix(&converted_ellpack_matrix);
    }
    fclose(resultsSer);
    fclose(resultsCsrOmp);
    fclose(resultsEllpackOmp);
    /*
    mat = getMatrix("test_matrices/matrix_files/prova_simm.mtx");
    
    matrix new_matrix;
    new_matrix.m = mat.m;
    new_matrix.n = mat.n;
    new_matrix.coeff = (double **)malloc(sizeof(double*)*mat.m);
    for(int i = 0; i < mat.m; i++){
        new_matrix.coeff[i] = (double *) calloc(mat.n, sizeof(double));     // allocazione della memoria per la riga i-esima
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
        new_matrix.coeff[mat.rows[i]][mat.cols[i]] = (double)(mat.values[i]);
    }
    stampaMatrice(new_matrix);
    */
}