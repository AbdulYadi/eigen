#include <stdio.h>
#include <math.h>
#include "eigen.h"

static void display(const matrix *m);

int main(int argc, char* argv[])
{
    matrix a, eigen_value, eigen_vector;
    int err;

    matrix_create(3, 3, &a);

    MATRIX(a, 0, 0) = 1;
    MATRIX(a, 0, 1) = 2;
    MATRIX(a, 0, 2) = 3;

    //second row
    MATRIX(a, 1, 0) = 4;
    MATRIX(a, 1, 1) = 5;
    MATRIX(a, 1, 2) = 6;

    //second row
    MATRIX(a, 2, 0) = 7;
    MATRIX(a, 2, 1) = 8;
    MATRIX(a, 2, 2) = 9;

    printf("Matrix A\n");
    display(&a);

    matrix_init(&eigen_value);
    matrix_init(&eigen_vector);
    err = eigen(&a, &eigen_value, &eigen_vector);
    printf("\n%s\n", eigen_err(err));
    if(err==EIGEN_OK)
    {
        printf("\nEigen values\n");
        display(&eigen_value);

        printf("\nEigen vector\n");
        display(&eigen_vector);
    }

    matrix_destroy(&a);
    matrix_destroy(&eigen_value);
    matrix_destroy(&eigen_vector);
    return 0;
}

static void display(const matrix *m)
{
    int r, c;
    for(r=0; r<m->rows; r++)
    {
        for(c=0; c<m->cols; c++)
            printf("%s%.04f", c>0 ? "\t" : "", MATRIXP(m, r, c));
        printf("\n");
    }
}
