#include <malloc.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "product/headers/product_csr.h"
#include "product/headers/product_ellpack_openmp.h"

#define NFILES 30
#define NREPETITIONS 2
#define FILEHEADER "nThreads, MatrixName, NRows, NCol, NZ, k, time, err_rel, GFLOPS\n"
#define MAX_N_THREADS 40

int main(){
    coo_matrix mat;
    csr_matrix converted_csr_matrix;
    ellpack_matrix converted_ellpack_matrix;
    matrix multivector, result_ser, result_par, multivector_T;
    double error_csr_omp, error_ellpack_omp, error_ellpack_cuda, error_csr_cuda;
    clock_t begin, end;
    double timeSumSer = 0.0, timeSumCsrOmp = 0.0, timeSumEllpackOmp = 0.0, timeSumCsrCuda = 0.0;
    int col_multivector[7] = {3,4,8,12,16,32,64};

    FILE *resultsSer = fopen("../measurements/results/serial.csv", "w+");
    FILE *resultsCsrOmp = fopen("../measurements/results/crsomp.csv", "w+");
    FILE *resultsEllpackOmp = fopen("../measurements/results/ellpackomp.csv", "w+");
    FILE *resultsCsrCuda = fopen("../measurements/results/csrcuda.csv", "w+");

    if(resultsSer == NULL || resultsCsrOmp == NULL || resultsEllpackOmp == NULL || resultsCsrCuda == NULL)
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
    fprintf(resultsCsrCuda,"%s",FILEHEADER);

    for(int i = 0; i < NFILES; i++){
        
        mat = get_matrix(matFiles[i]);
        printf("%s %d %d\n",matFiles[i],mat.m,mat.nz);
        converted_csr_matrix = convert_to_csr_from_coo(mat);
        converted_ellpack_matrix = convert_coo_to_ellpack(mat);
        for(int j = 0; j < 7; j++){
            printf("Multivettore %d\n",col_multivector[j]);
            multivector = generate_multivector(mat.n, col_multivector[j]);
            prepara_risultato(converted_csr_matrix.m, multivector.n,&result_ser);
            prepara_risultato(converted_csr_matrix.m, multivector.n, &result_par);
            timeSumSer = 0.0;
            timeSumCsrCuda = 0.0;
            error_csr_cuda = 0;
            error_ellpack_cuda = 0;
            
            for(int k = 0; k < NREPETITIONS; k++){
                timeSumSer += calcola_prodotto_seriale(converted_csr_matrix,multivector, &result_ser);
                
                prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                timeSumCsrCuda += calcola_prodotto_csr_cuda(converted_csr_matrix, multivector, &result_par);
                error_csr_cuda += check_result(result_ser,result_par);
                free_matrix(&result_par);
            }

            for(int t = 1; t < MAX_N_THREADS; t++){
                timeSumCsrOmp = 0.0;
                timeSumEllpackOmp = 0.0;  
                error_csr_omp = 0;  
                error_ellpack_omp = 0;
                for(int k = 0; k < NREPETITIONS; k++){
                    prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                    timeSumCsrOmp += calcola_prodotto_per_righe_csr_openmp(converted_csr_matrix,multivector, &result_par, t);
                    error_csr_omp += check_result(result_ser,result_par);
                    free_matrix(&result_par);

                    prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                    timeSumEllpackOmp += optimized_ellpack_product(converted_ellpack_matrix,multivector, &result_par, t);
                    error_ellpack_omp += check_result(result_ser,result_par);
                    free_matrix(&result_par);
                }   
                fprintf(resultsCsrOmp,"%d, %s, %d, %d, %d, %d, %f, %f, %f\n",t,matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumCsrOmp/NREPETITIONS,error_csr_omp/NREPETITIONS,(converted_csr_matrix.nz/pow(10,9))*(2*col_multivector[j]/(timeSumCsrOmp/NREPETITIONS)));
                fprintf(resultsEllpackOmp,"%d, %s, %d, %d, %d, %d, %f, %f, %f\n",t,matFiles[i],converted_ellpack_matrix.m, converted_ellpack_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumEllpackOmp/NREPETITIONS,error_ellpack_omp/NREPETITIONS,((2*col_multivector[j]/pow(10,9))*converted_csr_matrix.nz)/(timeSumEllpackOmp/NREPETITIONS));
            }

            fprintf(resultsSer,"%d, %s, %d, %d, %d, %d, %f, %f, %f\n",1, matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumSer/NREPETITIONS,0,(converted_csr_matrix.nz/pow(10,9))*(2*col_multivector[j]/(timeSumSer/NREPETITIONS)));
            fprintf(resultsCsrCuda,"%d, %s, %d, %d, %d, %d, %f, %f, %f\n",0,matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], timeSumCsrCuda/NREPETITIONS,error_csr_cuda/NREPETITIONS,((2*col_multivector[j]/pow(10,9))*converted_csr_matrix.nz)/(timeSumCsrCuda/NREPETITIONS));

            free_matrix(&result_ser);
        }
        free_matrix(&multivector);
        free_csr_matrix(&converted_csr_matrix);
        free_ellpack_matrix(&converted_ellpack_matrix);
    }

    fclose(resultsSer);
    fclose(resultsCsrOmp);
    fclose(resultsEllpackOmp);
    fclose(resultsCsrCuda);
    return 0;
}