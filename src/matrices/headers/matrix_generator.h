#ifndef _MATRIXGENERATORH_
#define _MATRIXGENERATORH_


#define MAX_RANDOM_VALUE 10

typedef struct{
    int m;              //Numero righe matrice
    int n;              //Numero colonne matrice
    int nz;             //Numero di non zeri della matrice
    double**      coeff; //Vettore dei coefficienti
} sparse_matrix;

typedef struct{
    int m;              //Numero righe multivettore
    int n;              //Numero colonne multivettore
    double*      coeff; //Vettore dei coefficienti
} matrix;

void stampa_matrice(matrix mat);

void stampa_matrice_su_file(char* filename, matrix mat);
void stampa_matrice_sparsa(sparse_matrix mat);

sparse_matrix generate_sparse_matrix(int m, int n, int max_nz);
matrix generate_multivector(int m, int n);

matrix genera_trasposta(matrix mat);
#endif