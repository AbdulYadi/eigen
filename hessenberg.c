#include <stdlib.h>
#include <math.h>
#include "matrix.h"

static void householder(matrix *h, int high, int m, matrix *ort, float scale);

void hessenberg(const matrix *m, matrix *h, matrix *eigen_vector)
{
    int low, high, _m, r, c;
    float scale, g;
    matrix ort;

    matrix_dup(m, h);
    matrix_identity_create_buffer(h->rows, eigen_vector);
    matrix_create(h->rows, 1, &ort);//single column matrix

    low = 0;
    high = m->rows - 1;

    for (_m = low + 1; _m <= high - 1; _m++)
    {
        scale = 0.0;
        for (r = _m; r <= high; r++)
        	scale += fabsf(MATRIXP(h, r, _m-1));

        if (scale != 0.0)
            householder(h, high, _m, &ort, scale);
    }

////Accumulate transformations (Algol's ortran)
    for (_m = high - 1; _m >= low + 1; _m--)
    {
        if(MATRIXP(h, _m, _m-1) != 0)
        {
            for (r = _m + 1; r <= high; r++)
                MATRIX(ort, r, 0) = MATRIXP(h, r, _m-1);

            for (c = _m; c <= high; c++) {
    			g = 0.0;
    			for (r = _m; r <= high; r++)
    				g += MATRIX(ort, r, 0) * MATRIXP(eigen_vector, r, c);

    			// Double division avoids possible underflow
    			g = (g / MATRIX(ort, _m, 0)) / MATRIXP(h, _m, _m-1);
    			for (r = _m; r <= high; r++)
                    MATRIXP(eigen_vector, r, c) += g * MATRIX(ort, r, 0);
    		}
        }
	}

    matrix_destroy(&ort);
}

static void householder(matrix *h, int high, int m, matrix *ort, float scale)
{
    float _h, g, f;
    int r, c;

    for (r = high; r >= m; r--)
    {
        MATRIXP(ort, r, 0) = MATRIXP(h, r, m-1) / scale;
		_h += MATRIXP(ort, r, 0) * MATRIXP(ort, r, 0);
    }

    g = sqrt(_h);
	if (MATRIXP(ort, m, 0) > 0)
		g = -g;

    _h -= MATRIXP(ort, m, 0) * g;
    MATRIXP(ort, m, 0) -= g;

    // Apply Householder similarity transformation
    // H = (I-u*u'/h)*H*(I-u*u')/h)
    for (c = m; c < h->cols; c++) {
		f = 0.0;
		for (r = high; r >= m; r--)
			f += MATRIXP(ort, r, 0) * MATRIXP(h, r, c);
		f = f / _h;
		for (r = m; r <= high; r++)
            MATRIXP(h, r, c) -= f * MATRIXP(ort, r, 0);
    }

    for (r = 0; r <= high; r++)
    {
		f = 0.0;
		for (c = high; c >= m; c--)
			f += MATRIXP(ort, c, 0) * MATRIXP(h, r, c);
		f = f / _h;
		for (c = m; c <= high; c++)
			MATRIXP(h, r, c) -= f * MATRIXP(ort, c, 0);
	}
	MATRIXP(ort, m, 0) = scale * MATRIXP(ort, m, 0);
    MATRIXP(h, m, m-1) = scale * g;
}
