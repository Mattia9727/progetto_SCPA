//
// Created by mattia971 on 3/28/23.
//
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define ROWS 10
#define COLS 10

int main() {
    int **sparse_matrix;    // doppio puntatore per la matrice
    int max_nz = 5;         // numero massimo di elementi non nulli per riga
    int row_nz;             // numero di elementi non nulli nella riga corrente
    int i, j, k;

    // Allocazione della memoria per il doppio puntatore
    sparse_matrix = (int **) malloc(sizeof(int *) * ROWS);
    if (sparse_matrix == NULL) {
        printf("Errore di allocazione della memoria.\n");
        return 1;
    }

    // Generazione casuale dei valori della matrice
    srand(time(NULL));
    for (i = 0; i < ROWS; i++) {
        row_nz = rand() % (max_nz + 1);     // numero casuale di elementi non nulli nella riga
        sparse_matrix[i] = (int *) calloc(COLS, sizeof(int));     // allocazione della memoria per la riga i-esima
        if (sparse_matrix[i] == NULL) {
            printf("Errore di allocazione della memoria.\n");
            return 1;
        }
        for (j = 0; j < row_nz; j++) {
            k = rand() % COLS;      // indice casuale di colonna per il valore non nullo
            sparse_matrix[i][k] = rand() % 10 + 1;     // valore casuale tra 1 e 10
        }
    }

    // Stampa della matrice
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++) {
            printf("%d ", sparse_matrix[i][j]);
        }
        printf("\n");
    }

    // Deallocazione della memoria
    for (i = 0; i < ROWS; i++) {
        free(sparse_matrix[i]);
    }
    free(sparse_matrix);

    return 0;
}

