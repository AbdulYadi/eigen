#include <stdlib.h>
#include <math.h>
#include "matrix.h"

static void householder(matrix *h, int high, int m, float *ort, float scale);

void hessenberg(const matrix *m, matrix *h, matrix *eigen_vector)
{
    int low, high, _m, r, c;
    float scale, g, *ort;

    matrix_dup(m, h);
    matrix_identity_create_buffer(h->rows, eigen_vector);

    low = 0;
    high = m->rows - 1;
    ort = (float*)malloc(m->rows * sizeof(float));//ort is column matrix

    for (_m = low + 1; _m <= high - 1; _m++)
    {
        scale = 0.0;
        for (r = _m; r <= high; r++)
        	scale += fabsf(MATRIXP(h, r, _m-1));

        if (scale != 0.0)
            householder(h, high, _m, ort, scale);
    }

////Accumulate transformations (Algol's ortran)
    for (_m = high - 1; _m >= low + 1; _m--)
    {
        if(MATRIXP(h, _m, _m-1) != 0)
        {
            for (r = _m + 1; r <= high; r++)
                ort[r] = MATRIXP(h, r, _m-1);

            for (c = _m; c <= high; c++) {
    			g = 0.0;
    			for (r = _m; r <= high; r++)
    				g += ort[r] * MATRIXP(eigen_vector, r, c);

    			// Double division avoids possible underflow
    			g = (g / ort[_m]) / MATRIXP(h, _m, _m-1);
    			for (r = _m; r <= high; r++)
                    MATRIXP(eigen_vector, r, c) += g * ort[r];
    		}
        }
	}
    free(ort);
}

static void householder(matrix *h, int high, int m, float *ort, float scale)
{
    float _h, g, f;
    int r, c;

    for (r = high; r >= m; r--)
    {
        ort[r] = MATRIXP(h, r, m-1) / scale;
		_h += ort[r] * ort[r];
    }

    g = sqrt(_h);
	if (ort[m] > 0)
		g = -g;

    _h -= ort[m] * g;
    ort[m] -= g;

    // Apply Householder similarity transformation
    // H = (I-u*u'/h)*H*(I-u*u')/h)
    for (c = m; c < h->cols; c++) {
		f = 0.0;
		for (r = high; r >= m; r--)
			f += ort[r] * MATRIXP(h, r, c);
		f = f / _h;
		for (r = m; r <= high; r++)
            MATRIXP(h, r, c) -= f * ort[r];
    }

    for (r = 0; r <= high; r++)
    {
		f = 0.0;
		for (c = high; c >= m; c--)
			f += ort[c] * MATRIXP(h, r, c);
		f = f / _h;
		for (c = m; c <= high; c++)
			MATRIXP(h, r, c) -= f * ort[c];
	}
	ort[m] = scale * ort[m];
    MATRIXP(h, m, m-1) = scale *g;
}
