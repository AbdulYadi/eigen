#include <stdlib.h>
#include <math.h>
#include <sys/param.h>
#include "matrix.h"

static void eigenvalue_outerloop(matrix *h, matrix *eigen_vector, int norm, matrix *eigen_value, matrix *e, float eps);
static void back_substitute(matrix *h, const matrix *eigen_value, const matrix *e, int norm, float eps);
static void cdiv(float xr, float xi, float yr, float yi, float* cdivr, float* cdivi);

void schur(matrix *h, matrix *eigen_vector, matrix *eigen_value)
{
    int low, high, row, col, k;
    float eps, norm, _z;
    matrix e;

    matrix_create(1, h->cols, eigen_value);

    low = 0;
    high = h->rows - 1;
    eps = pow(2.0, -52.0);
    matrix_create(h->rows, 1, &e);//single column matrix

    //store roots isolated by balance and compute matrix norm
    norm = 0.0;
    for (row = 0; row < h->rows; row++)
    {
    	if (row < low || row > high) {
            MATRIXP(eigen_value, 0, row) = MATRIXP(h, row, row);
            MATRIX(e, row, 0) = 0.0;
    	}
    	for (col = MAX(row - 1, 0); col < h->rows; col++)
    		norm += fabsf(MATRIXP(h, row, col));
    }

    //outer loop over eigenvalue index
    eigenvalue_outerloop(h, eigen_vector, norm, eigen_value, &e, eps);

    //backsubstitute to find vectors of upper triangular form
    if (norm == 0.0)
		return;
    back_substitute(h, eigen_value, &e, norm, eps);
    matrix_destroy(&e);

    //vectors of isolated roots
    for (row = 0; row < h->rows; row++)
    {
        if (row < low || row > high)
        {
            for (col = row; col < h->cols; col++)
                MATRIXP(eigen_vector, row, col) = MATRIXP(h, row, col);
        }
    }

    //back transformation to get eigenvectors of original matrix
    for (col = h->cols - 1; col >= low; col--)
    {
		for (row = low; row <= high; row++) {
			_z = 0.0;
			for (k = low; k <= MIN(col, high); k++)
                _z += MATRIXP(eigen_vector, row, k) * MATRIXP(h, k, col);
			MATRIXP(eigen_vector, row, col) = _z;
		}
    }
}

static void eigenvalue_outerloop(matrix *h, matrix *eigen_vector, int norm, matrix *eigen_value, matrix *e, float eps)
{
////assumption: h is squared matrix
    int n, m, k, iter, low, high, l, row, col;
    float _exshift, _w, _p, _q, _x, _y, _r, _s, _z;

    bool notlast;

    n = h->rows - 1;
    iter = 0;
    low = 0;
    high = h->rows - 1;
    _exshift = _w = _p = _q = _x = _y = _r = _s = _z = 0;

    _p = _q = 0.0;
	while (n >= low) {
		//look for single small sub-diagonal element
		l = n;
		while (l > low) {
			_s = fabsf(MATRIXP(h, l-1, l-1)) + fabsf(MATRIXP(h, l, l));
			if (_s == 0.0)
				_s = norm;
            if( fabsf(MATRIXP(h, l, l-1)) < eps * _s )
                break;
			l--;
		}

		//check for convergence
		if (l == n) {//one root found
            MATRIXP(h, n, n) += _exshift;
            MATRIXP(eigen_value, 0, n) = MATRIXP(h, n, n);
            MATRIXP(e, n, 0) = 0.0;
            n--;
            iter = 0;
		} else if (l == n - 1) {//two roots found
            _w = MATRIXP(h, n, n-1) * MATRIXP(h, n-1, n);
            _p = (MATRIXP(h, n-1, n-1) - MATRIXP(h, n, n)) / 2.0;
            _q = _p * _p + _w;
            _z = sqrt(fabsf(_q));
            MATRIXP(h, n, n) += _exshift;
            MATRIXP(h, n-1, n-1) += _exshift;
            _x = MATRIXP(h, n, n);

			//real pair
			if (_q >= 0) {
				if (_p >= 0)
					_z = _p + _z;
				else
					_z = _p - _z;
                MATRIXP(eigen_value, 0, n-1) =  _x + _z;
                MATRIXP(eigen_value, 0, n) = MATRIXP(eigen_value, 0, n-1);
				if (_z != 0.0)
					MATRIXP(eigen_value, 0, n) = _x - _w / _z;
				MATRIXP(e, n-1, 0) = 0.0;
				MATRIXP(e, n, 0) = 0.0;
                _x = MATRIXP(h, n, n-1);
				_s = fabsf(_x) + fabsf(_z);
				_p = _x / _s;
				_q = _z / _s;
				_r = sqrt(_p * _p + _q * _q);
				_p /= _r;
				_q /= _r;

				// Row modification
				for (col = n - 1; col < h->cols; col++) {
                    _z = MATRIXP(h, n-1, col);
                    MATRIXP(h, n-1, col) = _q * _z + _p * MATRIXP(h, n, col);
                    MATRIXP(h, n, col) = _q * MATRIXP(h, n, col) - _p * _z;
				}

				// Column modification
				for (row = 0; row <= n; row++) {
                    _z = MATRIXP(h, row, n-1);
                    MATRIXP(h, row, n-1) = _q * _z + _p * MATRIXP(h, row, n);
                    MATRIXP(h, row, n) = _q * MATRIXP(h, row, n) - _p * _z;
				}

				// Accumulate transformations
				for (row = low; row <= high; row++) {
                    _z = MATRIXP(eigen_vector, row, n-1);
                    MATRIXP(eigen_vector, row, n-1) = _q * _z + _p * MATRIXP(eigen_vector, row, n);
                    MATRIXP(eigen_vector, row, n) = _q * MATRIXP(eigen_vector, row, n) - _p * _z;
                }

			} else {// Complex pair
                MATRIXP(eigen_value, 0, n-1) = _x + _p;
                MATRIXP(eigen_value, 0, n) = _x + _p;
				MATRIXP(e, n-1, 0) = _z;
				MATRIXP(e, n, 0) = -_z;
			}
			n = n - 2;
			iter = 0;
		} else {// No convergence yet: no one or two root found
			// Form shift
            _x = MATRIXP(h, n, n);
			_y = 0.0;
			_w = 0.0;
			if (l < n) {
                _y = MATRIXP(h, n-1, n-1);
                _w = MATRIXP(h, n, n-1) * MATRIXP(h, n-1, n);
			}

			// Wilkinson's original ad hoc shift
			if (iter == 10) {
				_exshift += _x;
				for (row = low; row <= n; row++)
                    MATRIXP(h, row, row) -= _x;
                _s = fabsf(MATRIXP(h, n, n-1)) + fabsf(MATRIXP(h, n-1, n-2));
				_x = _y = 0.75 * _s;
				_w = -0.4375 * _s * _s;
			}

			// MATLAB's new ad hoc shift
			if (iter == 30) {
				_s = (_y - _x) / 2.0;
				_s = _s * _s + _w;
				if (_s > 0) {
					_s = sqrt(_s);
					if (_y < _x)
						_s = -_s;
					_s = _x - _w / ((_y - _x) / 2.0 + _s);
					for (row = low; row <= n; row++)
                        MATRIXP(h, row, row) -= _s;
					_exshift += _s;
					_x = _y = _w = 0.964;
				}
			}

			iter = iter + 1; // (Could check iteration count here.)

			// Look for two consecutive small sub-diagonal elements
			m = n - 2;
			while (m >= l) {
                _z = MATRIXP(h, m, m);
				_r = _x - _z;
				_s = _y - _z;

                _p = (_r * _s - _w) / MATRIXP(h, m+1, m) + MATRIXP(h, m, m+1);
                _q = MATRIXP(h, m+1, m+1) - _z - _r - _s;
                _r = MATRIXP(h, m+2, m+1);
				_s = fabsf(_p) + fabsf(_q) + fabsf(_r);
				_p /= _s;
				_q /= _s;
				_r /= _s;
				if (m == l)
					break;
                if (fabsf(MATRIXP(h, m, m-1)) * (fabsf(_q) + fabsf(_r))
                    < eps * (fabsf(_p)
                            * ( fabsf(MATRIXP(h, m-1, m-1))
                                    + fabsf(_z)
                                    + fabsf(MATRIXP(h, m+1, m+1))
                                )
                            )
                )
                    break;
				m--;
			}

			for (row = m + 2; row <= n; row++) {
                MATRIXP(h, row, row-2) = 0.0;
				if (row > m + 2)
                    MATRIXP(h, row, row-3) = 0.0;
			}

			// Double QR step involving rows l:n and columns m:n
			for (k = m; k <= n - 1; k++) {
				notlast = (k != n - 1);
				if (k != m) {
                    _p = MATRIXP(h, k, k-1);
					_q = MATRIXP(h, k+1, k-1);
                    _r = notlast ? MATRIXP(h, k + 2, k - 1) : 0.0;
					_x = fabsf(_p) + fabsf(_q) + fabsf(_r);
					if (_x != 0.0) {
						_p /= _x;
						_q /= _x;
						_r /= _x;
					}
				}
				if (_x == 0.0)
					break;
				_s = sqrt(_p * _p + _q * _q + _r * _r);
				if (_p < 0)
					_s = -_s;
				if (_s != 0) {
					if (k != m)
						MATRIXP(h, k, k - 1) = -_s * _x;
					else if (l != m)
                        MATRIXP(h, k, k - 1) = -MATRIXP(h, k, k - 1);
					_p += _s;
					_x = _p / _s;
					_y = _q / _s;
					_z = _r / _s;
					_q /= _p;
					_r /= _p;

					// Row modification
					for (col = k; col < h->cols; col++) {
                        _p = MATRIXP(h, k, col) + _q * MATRIXP(h, k+1, col);
						if (notlast) {
                            _p += _r * MATRIXP(h, k+2, col);
                            MATRIXP(h, k+2, col) -= _p * _z;
						}
                        MATRIXP(h, k, col) -= _p * _x;
                        MATRIXP(h, k+1, col) -= _p * _y;
					}

					// Column modification
					for (row = 0; row <= MIN(n, k + 3); row++) {
                        _p = _x * MATRIXP(h, row, k) + _y * MATRIXP(h, row, k+1);
						if (notlast) {
                            _p = _p + _z * MATRIXP(h, row, k+2);
                            MATRIXP(h, row, k+2) -= _p * _r;
						}
                        MATRIXP(h, row, k) -= _p;
                        MATRIXP(h, row, k+1) -= _p * _q;
					}

					// Accumulate transformations
					for (row = low; row <= high; row++) {
                        _p = _x * MATRIXP(eigen_vector, row, k) + _y * MATRIXP(eigen_vector, row, k+1);
						if (notlast) {
                            _p += _z * MATRIXP(eigen_vector, row, k+2);
                            MATRIXP(eigen_vector, row, k+2) -= _p * _r;
						}
                        MATRIXP(eigen_vector, row, k) -= _p;
                        MATRIXP(eigen_vector, row, k+1) -= _p * _q;
					}
				} // (s != 0)
			} // k loop
		} // check convergence
	} // while (n >= low)
}

static void back_substitute(matrix *h, const matrix *eigen_value, const matrix *e, int norm, float eps)
{
////assumption: h is squared matrix
    float _p, _q, _w, _x, _y, _t, _r, _s, _z, cdivr, cdivi, ra, sa, vr, vi;
    int n, l, row, col, row2;

    _p = _q = _w = _x = _y = _t = _r = _s = _z = 0;

    for (n = h->rows - 1; n >= 0; n--)
    {
		_p = MATRIXP(eigen_value, 0, n);
		_q = MATRIXP(e, n, 0);

		// Real vector
		if (_q == 0) {
			l = n;
            MATRIXP(h, n, n) = 1.0;
			for (row = n - 1; row >= 0; row--) {
                _w = MATRIXP(h, row, row) - _p;
				_r = 0.0;
				for (col = l; col <= n; col++)
                    _r += MATRIXP(h, row, col) * MATRIXP(h, col, n);
				if (MATRIXP(e, row, 0) < 0.0) {
					_z = _w;
					_s = _r;
				} else {
					l = row;
					if (MATRIXP(e, row, 0) == 0.0) {
						if (_w != 0.0)
                            MATRIXP(h, row, n) = -_r / _w;
						else
                            MATRIXP(h, row, n) = -_r / (eps * norm);
						// Solve real equations
					} else {
                        _x = MATRIXP(h, row, row+1);
                        _y = MATRIXP(h, row+1, row);
                        _q = (MATRIXP(eigen_value, 0, row) - _p) * (MATRIXP(eigen_value, 0, row) - _p) + MATRIXP(e, row, 0) * MATRIXP(e, row, 0);
						_t = (_x * _s - _z * _r) / _q;
                        MATRIXP(h, row, n) = _t;
						if (fabsf(_x) > fabsf(_z))
                            MATRIXP(h, row+1, n) = (-_r - _w * _t) / _x;
						else
                            MATRIXP(h, row+1, n) = (-_s - _y * _t) / _z;
					}

					// Overflow control
                    _t = fabsf(MATRIXP(h, row, n));
					if ((eps * _t) * _t > 1) {
						for (row2 = row; row2 <= n; row2++)
                            MATRIXP(h, row2, n) /= _t;
					}
				}
			}

			// Complex vector
		} else if (_q < 0) {
			int l = n - 1;

			// Last vector component imaginary so matrix is triangular
            if (fabsf( MATRIXP(h, n, n-1) ) > fabsf( MATRIXP(h, n-1, n) )) {
                MATRIXP(h, n-1, n-1) = _q / MATRIXP(h, n, n-1);
                MATRIXP(h, n-1, n) = -(MATRIXP(h, n, n) - _p) / MATRIXP(h, n, n-1);
			} else {
                cdiv(0.0, -MATRIXP(h, n-1, n), MATRIXP(h, n-1, n-1) - _p, _q, &cdivr, &cdivi);
                MATRIXP(h, n-1, n-1) = cdivr;
                MATRIXP(h, n-1, n) = cdivi;
			}
            MATRIXP(h, n, n-1) = 0.0;
            MATRIXP(h, n, n) = 1.0;
			for (row = n - 2; row >= 0; row--) {
				ra = 0.0;
				sa = 0.0;
				for (col = l; col <= n; col++) {
					ra = ra + MATRIXP(h, row, col) * MATRIXP(h, col, n-1);
					sa = sa + MATRIXP(h, row, col) * MATRIXP(h, col, n);
				}
				_w =  MATRIXP(h, row, row) - _p;

				if (MATRIXP(e, row, 0) < 0.0) {
					_z = _w;
					_r = ra;
					_s = sa;
				} else {
					l = row;
					if (MATRIXP(e, row, 0) == 0) {
                        cdiv(-ra, -sa, _w, _q, &cdivr, &cdivi);
                        MATRIXP(h, row, n-1) = cdivr;
                        MATRIXP(h, row, n) = cdivi;
					} else {
						// Solve complex equations
						_x = MATRIXP(h, row, row+1);
						_y = MATRIXP(h, row+1, row);
						vr = (MATRIXP(eigen_value, 0, row) - _p) * (MATRIXP(eigen_value, 0, row) - _p) + MATRIXP(e, row, 0) * MATRIXP(e, row, 0) - _q * _q;
						vi = (MATRIXP(eigen_value, 0, row) - _p) * 2.0 * _q;
						if (vr == 0.0 && vi == 0.0)
							vr = eps * norm * (fabsf(_w) + fabsf(_q) + fabsf(_x)
									+ fabsf(_y) + fabsf(_z));
						cdiv(_x * _r - _z * ra + _q * sa,
								_x * _s - _z * sa - _q * ra,
                                vr, vi, &cdivr, &cdivi);
                        MATRIXP(h, row, n-1) = cdivr;
						MATRIXP(h, row, n) = cdivi;
						if (fabsf(_x) > (fabsf(_z) + fabsf(_q))) {
                            MATRIXP(h, row+1, n-1) = (-ra - _w * MATRIXP(h, row, n - 1) + _q * MATRIXP(h, row, n)) / _x;
                            MATRIXP(h, row+1, n) = (-sa - _w * MATRIXP(h, row, n) - _q * MATRIXP(h, row, n-1)) / _x;
						} else {
							cdiv(-_r - _y * MATRIXP(h, row, n-1), -_s - _y * MATRIXP(h, row, n), _z, _q, &cdivr, &cdivi);
							MATRIXP(h, row+1, n - 1) = cdivr;
							MATRIXP(h, row+1, n ) = cdivi;
						}
					}

					// Overflow control
					_t = MAX(fabsf(MATRIXP(h, row, n-1)), fabsf(MATRIXP(h, row, n)));
					if ((eps * _t) * _t > 1) {
						for (row2 = row; row2 <= n; row2++) {
							MATRIXP(h, row2, n - 1) /= _t;
                            MATRIXP(h, row2, n) /= _t;
						}
					}
				}
			}
		}
	}//for (n = h->rows - 1; n >= 0; n--)
}

static void cdiv(float xr, float xi, float yr, float yi, float* cdivr, float* cdivi)
{
    float r, d;
    if (fabsf(yr) > fabsf(yi)) {
        r = yi / yr;
        d = yr + r * yi;
        *cdivr = (xr + r * xi) / d;
        *cdivi = (xi - r * xr) / d;
    } else {
        r = yr / yi;
        d = yi + r * yr;
        *cdivr = (r * xr + xi) / d;
        *cdivi = (r * xi - xr) / d;
    }
}
