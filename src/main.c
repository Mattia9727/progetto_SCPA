#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "product/headers/product_csr.h"
#include "product/headers/product_ellpack_openmp.h"

#define NFILES 1
#define NREPETITIONS 2
#define FILEHEADER "MatrixName, NRows, NCol, NZ, k, time, GFLOPS\n"

int main(){
    coo_matrix mat;
    csr_matrix converted_csr_matrix;
    ellpack_matrix converted_ellpack_matrix;
    matrix multivector, result1, result2, multivector_T;
    clock_t begin, end;
    double timeSumSer = 0.0, timeSumCsrOmp = 0.0, timeSumEllpackOmp = 0.0;
    int col_multivector[7] = {3,4,8,12,16,32,64};

    FILE *resultsSer = fopen("../measurements/results/serial.csv", "w+");
    FILE *resultsCsrOmp = fopen("../measurements/results/crsomp.csv", "w+");
    FILE *resultsEllpackOmp = fopen("../measurements/results/ellpackomp.csv", "w+");

    if(resultsSer == NULL || resultsCsrOmp == NULL || resultsEllpackOmp == NULL)
        printf("errore apertura file\n");

    char* matFiles[] = {
        "../measurements/matrix_files/cage4.mtx",
        "../measurements/matrix_files/olm1000.mtx",
        "../measurements/matrix_files/west2021.mtx",
        "../measurements/matrix_files/mhda416.mtx",
        "../measurements/matrix_files/adder_dcop_32.mtx",
        "../measurements/matrix_files/mcfe.mtx",
        "../measurements/matrix_files/rdist2.mtx",
        "../measurements/matrix_files/cavity10.mtx",
        "../measurements/matrix_files/mhd4800a.mtx",
        "../measurements/matrix_files/bcsstk17.mtx",
        "../measurements/matrix_files/raefsky2.mtx",
        "../measurements/matrix_files/thermal1.mtx",
        "../measurements/matrix_files/af23560.mtx",
        "../measurements/matrix_files/thermomech_TK.mtx",
        "../measurements/matrix_files/olafu.mtx",
        "../measurements/matrix_files/FEM_3D_thermal1.mtx",
        "../measurements/matrix_files/lung2.mtx",
        "../measurements/matrix_files/dc1.mtx",
        "../measurements/matrix_files/amazon0302.mtx",
        "../measurements/matrix_files/roadNet-PA.mtx",
        "../measurements/matrix_files/cop20k_A.mtx",
        "../measurements/matrix_files/mac_econ_fwd500.mtx",
        "../measurements/matrix_files/cant.mtx",
        "../measurements/matrix_files/webbase-1M.mtx",
        "../measurements/matrix_files/thermal2.mtx",
        "../measurements/matrix_files/nlpkkt80.mtx",
        "../measurements/matrix_files/PR02R.mtx",
        "../measurements/matrix_files/af_1_k101.mtx",
        "../measurements/matrix_files/ML_Laplace.mtx",
        "../measurements/matrix_files/Cube_Coup_dt0.mtx",
    }; 

    fprintf(resultsSer,"%s",FILEHEADER);
    fprintf(resultsCsrOmp,"%s",FILEHEADER);
    fprintf(resultsEllpackOmp,"%s",FILEHEADER);

    for(int i = 0; i < NFILES; i++){
        printf("%s\n",matFiles[i]);
        mat = get_matrix(matFiles[i]);
        converted_csr_matrix = convert_to_csr_from_coo(mat);
        converted_ellpack_matrix = convert_coo_to_ellpack(mat);
        for(int j = 0; j < 7; j++){
            multivector = generate_multivector(mat.n, col_multivector[j]);
            //multivector_T = genera_trasposta(multivector);
            //printf("Trasposto\n");
            prepara_risultato(converted_csr_matrix.m, multivector.n,&result1);
            prepara_risultato(converted_csr_matrix.m, multivector.n, &result2);
            timeSumSer = 0.0;
            timeSumCsrOmp = 0.0;
            timeSumEllpackOmp = 0.0;
            for(int k = 0; k < NREPETITIONS; k++){
                printf("seriale\n");
                timeSumSer += calcola_prodotto_seriale(converted_csr_matrix,multivector, &result1);
                printf("csr\n");
                timeSumCsrOmp += calcola_prodotto_per_righe_csr_openmp(converted_csr_matrix,multivector, &result2);
                //free_matrix(&result1);
                //prepara_risultato(converted_csr_matrix.m, multivector.n,&result1);
                printf("ellpack\n");
                timeSumEllpackOmp += optimized_ellpack_product(converted_ellpack_matrix,multivector, &result1);
                //printf("%d\n",k);
                check_result(result1,result2);
                calcola_prodotto_csr_cuda(converted_csr_matrix, multivector, &result1);
            }

            fprintf(resultsSer,"%s, %d, %d, %d, %d, %f, %f\n",matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumSer/NREPETITIONS,(converted_csr_matrix.nz/pow(10,9))*(2*col_multivector[j]/(timeSumSer/NREPETITIONS)));
            fprintf(resultsCsrOmp,"%s, %d, %d, %d, %d, %f, %f\n",matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumCsrOmp/NREPETITIONS,(converted_csr_matrix.nz/pow(10,9))*(2*col_multivector[j]/(timeSumCsrOmp/NREPETITIONS)));
            fprintf(resultsEllpackOmp,"%s, %d, %d, %d, %d, %f, %f\n",matFiles[i],converted_ellpack_matrix.m, converted_ellpack_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumEllpackOmp/NREPETITIONS,((2*col_multivector[j]/pow(10,9))*converted_csr_matrix.nz)/(timeSumEllpackOmp/NREPETITIONS));

            free_matrix(&result1);
            free_matrix(&result2);
        }
        free_matrix(&multivector);
        free_csr_matrix(&converted_csr_matrix);
        free_ellpack_matrix(&converted_ellpack_matrix);
    }

    fclose(resultsSer);
    fclose(resultsCsrOmp);
    fclose(resultsEllpackOmp);
    return 0;
}