typedef struct{
    int m;              //Numero righe matrice
    int n;              //Numero colonne matrice
    int nz;             //Numero di non zeri della matrice
    float**      coeff; //Vettore dei coefficienti
} sparse_matrix;

typedef struct{
    int m;              //Numero righe multivettore
    int n;              //Numero colonne multivettore
    float**     coeff;  //Vettore dei coefficienti
} multivector;