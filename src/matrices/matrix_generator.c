#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "headers/matrix_generator.h"

void stampa_matrice(matrix mat){
    int m = mat.m;
    int n = mat.n;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", mat.coeff[i * mat.n + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void stampa_matrice_su_file(char* filename, matrix mat){
    int m = mat.m;
    int n = mat.n;
    FILE *f = fopen(filename, "a");
    printf("file aperto\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(f,"%f ", mat.coeff[i * mat.n + j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

void stampa_matrice_sparsa(sparse_matrix mat){
    int m = mat.m;
    int n = mat.n;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (mat.coeff[i][j]==0.0) printf("\033[0;33m"); //Set the text to the color red
            printf("%f ", mat.coeff[i][j]);
            if (mat.coeff[i][j]==0.0) printf("\033[0m"); //Resets the text to default color

        }
        printf("\n");
    }
}

sparse_matrix generate_sparse_matrix(int m, int n, int max_nz) {
    sparse_matrix new_matrix;
    new_matrix.m = m;
    new_matrix.n = n;
    new_matrix.nz = 0;

    int row_nz;             // numero di elementi non nulli nella riga corrente
    int i, j, k;

    // Allocazione della memoria per il doppio puntatore
    new_matrix.coeff = (double **) malloc(sizeof(double *) * m);
    if (new_matrix.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }

    // Generazione casuale dei valori della matrice
    srand(time(NULL));
    for (i = 0; i < m; i++) {
        row_nz = rand() % (max_nz + 1);     // numero casuale di elementi non nulli nella riga
        new_matrix.coeff[i] = (double *) calloc(n, sizeof(double));     // allocazione della memoria per la riga i-esima
        if (new_matrix.coeff[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            exit(0);
        }
        for (j = 0; j < row_nz; j++) {
            k = rand() % n;      // indice casuale di colonna per il valore non nullo
            new_matrix.coeff[i][k] = (double)rand()/(double)(RAND_MAX/MAX_RANDOM_VALUE);;     // valore casuale tra 1 e 10
        }
    }
    // Stampa della matrice
    // stampa_matrice_sparsa(new_matrix);

    for(int i=0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (new_matrix.coeff[i][j] !=0) new_matrix.nz++;
        }
    }

    return new_matrix;
}

matrix generate_multivector(int m, int n) {
    matrix new_multivector;
    new_multivector.m = m;
    new_multivector.n = n;
    new_multivector.coeff = (double *)calloc(m*n,sizeof(double));
    if (new_multivector.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }

    int i, j, k;
    
    // Generazione casuale dei valori della matrice
    srand(time(NULL));
    for (i = 0; i < m; i++) {

        for (j = 0; j < n; j++) {
            new_multivector.coeff[i*n + j] = (double)rand()/(double)(RAND_MAX/MAX_RANDOM_VALUE);;     // valore casuale tra 1 e 10
                   
        }
    }
    // Stampa della matrice
    //stampa_matrice(new_multivector);

    return new_multivector;
}

matrix genera_trasposta(matrix mat){
    matrix mTrasposta;
    mTrasposta.m = mat.n;
    mTrasposta.n = mat.m;
    mTrasposta.coeff = (double*)malloc(sizeof(double)*mat.m*mat.n);
    if (mTrasposta.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }

    for(int i = 0; i < mat.m; i++){
        for(int j = 0; j < mat.n; j++){
            mTrasposta.coeff[j*mat.m + i] =  mat.coeff[i * mat.n + j];
        }
    }

    return mTrasposta;

}