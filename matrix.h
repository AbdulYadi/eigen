#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdbool.h>

typedef struct matrix {
    unsigned int rows;
    unsigned int cols;
    float *data;
} matrix;

#define MATRIX(M,ROW,COL) M.data[(ROW) * M.cols + (COL)]
#define MATRIXP(PM,ROW,COL) (PM)->data[(ROW) * (PM)->cols + (COL)]

extern void matrix_create(unsigned int rows, unsigned int cols, matrix* m);
extern bool matrix_is_created(const matrix* m);
extern bool matrix_is_same_dimension(const matrix* m1, const matrix* m2);
extern bool matrix_is_same_rows(const matrix* m1, const matrix* m2);
extern void matrix_init(matrix* m);
extern void matrix_destroy(matrix* m);
extern void matrix_mul(const matrix* m1, const matrix* m2, matrix* m, bool create_buffer);
extern void matrix_scalar_mul(matrix* m, float scalar);
extern void matrix_sub(const matrix* m1, const matrix* m2, matrix* m, bool create_buffer);
extern void matrix_dup(const matrix* m, matrix *dup);
extern void matrix_identity(matrix* m);
extern void matrix_identity_create_buffer(int rows, matrix* m);
extern bool matrix_is_identity(const matrix* m, float tau);
extern bool matrix_is_squared(const matrix *m);
extern bool matrix_is_row(const matrix *m);
extern bool matrix_is_column(const matrix *m);
extern bool matrix_is_symmetric(const matrix* m);
extern void matrix_copy(const matrix* m, matrix *copy);
extern void matrix_transpose(const matrix* m, matrix *out, bool create_buffer);
extern void matrix_zero(matrix* m);
extern void matrix_swap_row(matrix* m, int r1, int r2);
extern void matrix_diagonal_create_buffer(int rows, float value, matrix* m);
extern void matrix_diagonal(matrix* m, float value);
extern void matrix_diagonal_sub(matrix* m, float value);
extern void matrix_diagonal_add(matrix* m, float value);
extern void matrix_diagonal_extract(const matrix* m, matrix* row, bool create_buffer);
extern void matrix_diagonal_expand(const matrix* row, matrix* square, bool create_buffer);

#endif //_MATRIX_H_
