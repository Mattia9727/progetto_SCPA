#include <malloc.h>
#include <stdio.h>
#include <stdlib.h>
#include "product/headers/product_csr.h"
#include "product/headers/product_ellpack.h"
#include <time.h>
#include <math.h>
#include "matrices/format/headers/ellpack.h"

#define NFILES 30
#define NREPETITIONS 10
#define FILEHEADER "nThreads, rep, MatrixName, NRows, NCol, NZ, k, hackSize, time, err_rel, speedup, bandwidth, GFLOPS\n"
#define MAX_N_THREADS 40
#define STARTFILE 0

int main(){
    coo_matrix mat;
    csr_matrix converted_csr_matrix;
    ellpack_matrix converted_ellpack_matrix;
    h_ellpack_matrix_bis converted_h_ellpack_matrix_bis;
    matrix multivector, result_par, result_ser,multivector_T;
    performance perf;

    clock_t begin, end;
    double timeSumSer = 0.0, timeSumCsrOmp = 0.0, timeSumEllpackOmp = 0.0, timeSumEllpackCuda = 0.0, timeSumCsrCuda = 0.0;
    double errorCsrOmp = 0.0, errorEllpackOmp = 0.0, errorEllpackCuda = 0.0, errorCsrCuda = 0.0;
    double speedupCsrOmp = 0.0, speedupCsrCuda = 0.0, speedupEllpackOmp = 0.0, speedupEllpackCuda = 0.0;
    double bandwidthCsrOmp = 0.0, bandwidthCsrCuda = 0.0, bandwidthEllpackOmp = 0.0, bandwidthEllpackCuda = 0.0;
    int col_multivector[7] = {3,4,8,12,16,32,64};

    FILE *resultsSer = fopen("../measurements/results/serial.csv", "w+");
    FILE *resultsCsrOmp = fopen("../measurements/results/csromp.csv", "w+");
    FILE *resultsEllpackOmp = fopen("../measurements/results/ellpackomp.csv", "w+");
    FILE *resultsCsrCuda = fopen("../measurements/results/csrcuda.csv", "w+");
    FILE *resultsEllpackCuda = fopen("../measurements/results/ellpackcuda.csv", "w+");
    FILE *debug = fopen("../measurements/results/debug.csv", "w+");

    if(resultsSer == NULL || resultsCsrOmp == NULL || resultsEllpackOmp == NULL ||resultsCsrCuda == NULL || resultsEllpackCuda == NULL)
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
    fprintf(resultsEllpackCuda,"%s",FILEHEADER);

    float memory_clock_rate = 7.001; //GHz
    int memory_bus_width = 256; //bit
    float peak_theoretical_memory_bandwidth = (memory_clock_rate * (10^9) * (256/8) * 2) / (10^9);
    printf("Banda teorica GPU: %f GBps\n\n", peak_theoretical_memory_bandwidth);
    printf("Banda teorica CPU: 204.8 GBps\n\n");

    for(int i = STARTFILE; i < NFILES; i++){
        printf("numero %d, matrice %s\n",i,matFiles[i]);
        mat = get_matrix(matFiles[i]);
        converted_csr_matrix = convert_to_csr_from_coo(mat);
        converted_ellpack_matrix = convert_coo_to_ellpack(mat);
        
        for(int j = 0; j < 7; j++){
            multivector = generate_multivector(mat.n, col_multivector[j]);
            prepara_risultato(converted_csr_matrix.m, multivector.n, &result_ser);
            prepara_risultato(converted_csr_matrix.m, multivector.n, &result_par);
            

            for(int k = 0; k < NREPETITIONS; k++){
                timeSumSer = calcola_prodotto_seriale(converted_csr_matrix,multivector, &result_ser);
            
                prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                perf = calcola_prodotto_csr_cuda(converted_csr_matrix, multivector, &result_par);
                timeSumCsrCuda = perf.time;
                bandwidthCsrCuda = perf.bandwidth;
                errorCsrCuda = check_result(result_ser,result_par);
                speedupCsrCuda = timeSumSer/timeSumCsrCuda;
                free_matrix(&result_par);
                
                fprintf(resultsSer,"%d, %d, %s, %d, %d, %d, %d, %d, %f, %f, 1, %f, %f\n",1, k, matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], 0,timeSumSer,0,0,(converted_csr_matrix.nz/pow(10,9))*(2*col_multivector[j]/(timeSumSer)));
                fprintf(resultsCsrCuda,"%d, %d, %s, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f\n",0,k,matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], 0,timeSumCsrCuda,errorCsrCuda,speedupCsrCuda,bandwidthCsrCuda,((2*col_multivector[j]/pow(10,9))*converted_csr_matrix.nz)/(timeSumCsrCuda));
             
            }
            
            for(int hackOffset = 32; hackOffset <= 64; hackOffset+=32){
                converted_h_ellpack_matrix_bis = convert_coo_to_h_ellpack_bis(mat,hackOffset);
                for(int k = 0; k < NREPETITIONS; k++){
                    prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                    perf = optimized_cuda_h_ellpack_product_in_bis(converted_h_ellpack_matrix_bis, multivector, &result_par);
                    timeSumEllpackCuda = perf.time;
                    bandwidthEllpackCuda = perf.bandwidth;
                    errorEllpackCuda = check_result(result_ser,result_par);
                    speedupEllpackCuda = timeSumSer/timeSumEllpackCuda;
                    free_matrix(&result_par);
                    fprintf(resultsEllpackCuda,"%d, %d, %s, %d, %d, %d, %d, %d, %f, %f, %f, %f, %f\n",0,k,matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], hackOffset, timeSumEllpackCuda,errorEllpackCuda,speedupEllpackCuda,bandwidthEllpackCuda,((2*col_multivector[j]/pow(10,9))*converted_csr_matrix.nz)/(timeSumEllpackCuda));
                   
                }  
                free_h_ellpack_matrix(&converted_h_ellpack_matrix_bis);  
            }
            
            for(int t = 1; t < MAX_N_THREADS; t++){
                for(int k = 0; k < NREPETITIONS; k++){
                    prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                    timeSumCsrOmp = calcola_prodotto_per_righe_csr_openmp(converted_csr_matrix,multivector, &result_par, t);
                    errorCsrOmp = check_result(result_ser,result_par);
                    speedupCsrOmp = timeSumSer/timeSumCsrOmp;
                    free_matrix(&result_par);
                    
                    prepara_risultato(converted_csr_matrix.m, multivector.n,&result_par);
                    
                    timeSumEllpackOmp = optimized_omp_ellpack_product(converted_ellpack_matrix,multivector, &result_par, t);
                   
                    errorEllpackOmp = check_result(result_ser,result_par);
                    speedupEllpackOmp = timeSumSer/timeSumEllpackOmp;
                    free_matrix(&result_par);
                    fprintf(resultsCsrOmp,"%d, %d, %s, %d, %d, %d, %d, %d, %f, %f, %f, %f\n",t,k,matFiles[i],converted_csr_matrix.m, converted_csr_matrix.n, converted_csr_matrix.nz, col_multivector[j], 0, timeSumCsrOmp,errorCsrOmp,speedupCsrOmp,0,(converted_csr_matrix.nz/pow(10,9))*(2*col_multivector[j]/(timeSumCsrOmp)));
                    fprintf(resultsEllpackOmp,"%d, %d, %s, %d, %d, %d, %d, %d, %f, %f, %f, %f\n",t,k,matFiles[i],converted_ellpack_matrix.m, converted_ellpack_matrix.n, converted_csr_matrix.nz, col_multivector[j], 0, timeSumEllpackOmp,errorEllpackOmp,speedupEllpackOmp,0,((2*col_multivector[j]/pow(10,9))*converted_csr_matrix.nz)/(timeSumEllpackOmp));

                }
            }
            
            free_matrix(&result_ser);
            
        }
        
        free_coo_matrix(&mat);
        free_matrix(&multivector);        
    }

    fclose(resultsSer);
    fclose(resultsCsrOmp);
    fclose(resultsEllpackOmp);
    return 1;
}