#include "eigen.h"
#include "hessenberg.h"
#include "schur.h"

int eigen(const matrix* m, matrix* eigen_value, matrix *eigen_vector)
{
    matrix h;
    int r, c;
    float denom;

    if(!matrix_is_squared(m))
        return EIGEN_INPUT_NOT_SQUARED;

    matrix_init(&h);

    hessenberg(m, &h, eigen_vector);//eigen_vector is allocated in hessenberg
    schur(&h, eigen_vector, eigen_value);//eigen_value is allocated in schur

    //divide each column by last row then set last row to 1
    for(c=0; c<eigen_vector->cols; c++)
    {
        denom = MATRIXP(eigen_vector, eigen_vector->rows-1, c);
        if(denom!=1)
        {
            for(r=0; r<eigen_vector->rows-1; r++)
                MATRIXP(eigen_vector, r, c) /= denom;
            MATRIXP(eigen_vector, eigen_vector->rows-1, c) = 1;
        }
    }

    matrix_destroy(&h);
    return EIGEN_OK;
}

const char* eigen_err(int err)
{
    switch(err)
    {
        case EIGEN_OK:
            return "Success";
        case EIGEN_INPUT_NOT_SQUARED:
            return "Input matrix is not squared";
        default:
            break;
    }
    return "Unknown";
}
