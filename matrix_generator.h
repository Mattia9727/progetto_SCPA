#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX_RANDOM_VALUE 10

typedef struct{
    int m;              //Numero righe matrice
    int n;              //Numero colonne matrice
    int nz;             //Numero di non zeri della matrice
    float**      coeff; //Vettore dei coefficienti
} sparse_matrix;

typedef struct{
    int m;              //Numero righe multivettore
    int n;              //Numero colonne multivettore
    float**      coeff; //Vettore dei coefficienti
} multivector;



sparse_matrix GenerateSparseMatrix(int m, int n) {
    sparse_matrix new_matrix;
    new_matrix.m = m;
    new_matrix.n = n;
    new_matrix.nz = 0;
    new_matrix.coeff = (float **)malloc(sizeof(float*)*m);

    int max_nz = 5;         // numero massimo di elementi non nulli per riga
    int row_nz;             // numero di elementi non nulli nella riga corrente
    int i, j, k;

    // Allocazione della memoria per il doppio puntatore
    new_matrix.coeff = (float **) malloc(sizeof(float *) * m);
    if (new_matrix.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }

    // Generazione casuale dei valori della matrice
    srand(time(NULL));
    for (i = 0; i < m; i++) {
        row_nz = rand() % (max_nz + 1);     // numero casuale di elementi non nulli nella riga
        new_matrix.coeff[i] = (float *) calloc(n, sizeof(float));     // allocazione della memoria per la riga i-esima
        if (new_matrix.coeff[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            exit(0);
        }
        for (j = 0; j < row_nz; j++) {
            k = rand() % n;      // indice casuale di colonna per il valore non nullo
            new_matrix.coeff[i][k] = (float)rand()/(float)(RAND_MAX/MAX_RANDOM_VALUE);;     // valore casuale tra 1 e 10
        }
    }

    // Stampa della matrice
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", new_matrix.coeff[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");

    for(int i=0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (new_matrix.coeff[i][j] !=0) new_matrix.nz++;
        }
    }

    return new_matrix;
}

multivector GenerateMultivector(int m, int n) {
    multivector new_multivector;
    new_multivector.m = m;
    new_multivector.n = n;
    new_multivector.coeff = (float **)malloc(sizeof(float*)*m);

    int max_nz = 5;         // numero massimo di elementi non nulli per riga
    int i, j, k;

    // Allocazione della memoria per il doppio puntatore
    new_multivector.coeff = (float **) malloc(sizeof(float *) * m);
    if (new_multivector.coeff == NULL) {
        printf("Errore di allocazione della memoria.\n");
        exit(0);
    }

    // Generazione casuale dei valori della matrice
    srand(time(NULL));
    for (i = 0; i < m; i++) {
        new_multivector.coeff[i] = (float *) calloc(n, sizeof(float));     // allocazione della memoria per la riga i-esima
        if (new_multivector.coeff[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            exit(0);
        }
        for (j = 0; j < n; j++) {
            new_multivector.coeff[i][j] = (float)rand()/(float)(RAND_MAX/MAX_RANDOM_VALUE);;     // valore casuale tra 1 e 10
        }
    }

    // Stampa della matrice
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", new_multivector.coeff[i][j]);
        }
        printf("\n");
    }

    return new_multivector;
}