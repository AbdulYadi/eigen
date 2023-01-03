#ifndef _EIGEN_H_
#define _EIGEN_H_

#define EIGEN_OK 0
#define EIGEN_INPUT_NOT_SQUARED 1

#include "matrix.h"

extern int eigen(const matrix* m, matrix* eigen_value, matrix *eigen_vector);
extern const char* eigen_err(int err);

#endif //_EIGEN_H_
