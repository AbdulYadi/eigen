#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "matrix.h"

void matrix_create(unsigned int rows, unsigned int cols, matrix* m)
{
    m->rows = rows;
    m->cols = cols;
    m->data = (float*)malloc(rows * cols * sizeof(float));
}

bool matrix_is_created(const matrix* m)
{
    return m!=NULL && m->data!=NULL;
}

bool matrix_is_same_dimension(const matrix* m1, const matrix* m2)
{
    return m1->rows==m2->rows && m1->cols==m2->cols;
}

bool matrix_is_same_rows(const matrix* m1, const matrix* m2)
{
    return m1->rows==m2->rows;
}

void matrix_init(matrix* m)
{
    m->rows = m->cols = 0;
    m->data = NULL;
}

void matrix_destroy(matrix* m)
{
    if(m->data!=NULL) {
        free(m->data);
        m->data = NULL;
    }
    return;
}

void matrix_mul(const matrix* m1, const matrix* m2, matrix* m, bool create_buffer)
{
    int r, c, k;
    if(create_buffer)
        matrix_create(m1->rows, m2->cols, m);
    for(r=0; r<m1->rows; r++)
    {
        for(c=0; c<m2->cols; c++)
        {
            MATRIXP(m, r, c) = 0;
            for(k=0; k<m1->cols; k++)
                MATRIXP(m, r, c) += MATRIXP(m1, r, k) * MATRIXP(m2, k, c);
        }
    }
}

void matrix_scalar_mul(matrix* m, float scalar)
{
    int r, c;
    for(r=0; r<m->rows; r++)
    {
        for(c=0; c<m->cols; c++)
            MATRIXP(m, r, c) *= scalar;
    }
}

void matrix_sub(const matrix* m1, const matrix* m2, matrix* m, bool create_buffer)
{
    int r, c, rmin, cmin;

    rmin = m1->rows < m2->rows ? m1->rows : m2->rows;
    cmin = m1->cols < m2->cols ? m1->cols : m2->cols;

    if(!create_buffer)
    {
        rmin = m->rows < rmin ? m->rows : rmin;
        cmin = m->cols < cmin ? m->cols : cmin;
    }
    else
        matrix_create(rmin, cmin, m);

    for(r=0; r<rmin; r++)
    {
        for(c=0; c<cmin; c++)
            MATRIXP(m, r, c) = MATRIXP(m1, r, c) - MATRIXP(m2, r, c);
    }
}

void matrix_dup(const matrix* m, matrix *dup)
{
    matrix_create(m->rows, m->cols, dup);
    matrix_copy(m, dup);
}

void matrix_copy(const matrix* m, matrix *copy)
{
    memcpy(copy->data, m->data, m->rows * m->cols * sizeof(float));
}

void matrix_identity(matrix* m)
{
    matrix_diagonal(m, 1);
}

void matrix_identity_create_buffer(int rows, matrix* m)
{
    matrix_create(rows, rows, m);
    matrix_identity(m);
}

bool matrix_is_identity(const matrix* m, float tau)
{
    int r, c;
    float x;
    if(m->rows!=m->cols)
        return false;
    for(r=0; r<m->rows; r++)
    {
        for(c=0; c<m->cols; c++)
        {
            x = fabsf( MATRIXP(m, r, c) );
            if(r==c && fabsf(1 - x) > tau)
                return false;
            else if(r!=c && fabsf(0 - x) > tau)
                return false;
        }
    }
    return true;
}

bool matrix_is_squared(const matrix *m)
{
    return m->rows==m->cols;
}

bool matrix_is_row(const matrix *m)
{
    return m->rows==1 && m->cols>0;
}

bool matrix_is_column(const matrix *m)
{
    return m->rows>0 && m->cols==1;
}

bool matrix_is_symmetric(const matrix* m)
{
    int r, c;

    if(m->rows!=m->cols)
        return false;

    for(r=1; r<m->rows; r++)
    {
        for(c=0; c<r; c++)
        {
            if(MATRIXP(m, r, c)!=MATRIXP(m, c, r))
                return false;
        }
    }

    return true;
}

void matrix_transpose(const matrix* m, matrix *out, bool create_buffer)
{
    int r, c;
    if(create_buffer)
        matrix_create(m->cols, m->rows, out);
    for(r=0; r<m->rows; r++)
    {
        for(c=0; c<m->cols; c++)
            MATRIXP(out, c, r) = MATRIXP(m, r, c);
    }
}

void matrix_zero(matrix* m)
{
    memset(m->data, 0, m->rows * m->cols * sizeof(float));
}

void matrix_swap_row(matrix* m, int r1, int r2)
{
    float *swap;
    
    if(r1<0 || r1>m->rows-1 || r2<0 || r2>m->rows-1 || r1==r2)
        return;

    swap = (float*)malloc(m->cols * sizeof(float));
    memcpy(swap, &m->data[r1 * m->cols], m->cols * sizeof(float));
    memcpy(&m->data[r1 * m->cols], &m->data[r2 * m->cols], m->cols * sizeof(float));
    memcpy(&m->data[r2 * m->cols], swap, m->cols * sizeof(float));
    free(swap);
}

void matrix_diagonal_create_buffer(int rows, float value, matrix* m)
{
    matrix_create(rows, rows, m);
    matrix_diagonal(m, value);
}

void matrix_diagonal(matrix* m, float value)
{
    int i, min;
    matrix_zero(m);
    min = m->rows < m->cols ? m->rows : m->cols;
    for(i=0; i<min; i++)
        MATRIXP(m, i, i) = value;
}

void matrix_diagonal_sub(matrix* m, float value)
{
    int i, min;
    min = m->rows < m->cols ? m->rows : m->cols;
    for(i=0; i<min; i++)
        MATRIXP(m, i, i) -= value;
}

void matrix_diagonal_add(matrix* m, float value)
{
    int i, min;
    min = m->rows < m->cols ? m->rows : m->cols;
    for(i=0; i<min; i++)
        MATRIXP(m, i, i) += value;
}

void matrix_diagonal_extract(const matrix* m, matrix* row, bool create_buffer)
{
    int i, min;
    min = m->rows < m->cols ? m->rows : m->cols;
    if(!create_buffer)
        min = row->cols < min ? row->cols : min;
    else
        matrix_create(1, min, row);
    for(i=0; i<min; i++)
        row->data[i] = MATRIXP(m, i, i);
}

void matrix_diagonal_expand(const matrix* row, matrix* square, bool create_buffer)
{
    int i, min;
    min = row->cols;
    if(!create_buffer)
    {
        min = square->rows < min ? square->rows : min;
        min = square->cols < min ? square->cols : min;
    }
    else
        matrix_create(min, min, square);
    matrix_zero(square);
    for(i=0; i<min; i++)
        MATRIXP(square, i, i) = row->data[i];
}
