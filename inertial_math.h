//
// Created by tdhuang on 2022/4/12.
//

#ifndef GRIZZ_INS_INERTIAL_MATH_H
#define GRIZZ_INS_INERTIAL_MATH_H

#include <stddef.h>
#include <malloc.h>

enum{ERR_MAT_DIM=-1,ERR_NORM_ZERO=-2};

#define MID(C, i, j) ((i)*(C)+(j))

#define SQR(d) ((d)*(d))

/// zero-judge
/// \param d
/// \return
extern int is_zero(double d);

/// matrix zero-judge
/// \param d
/// \param n
/// \return
extern int all_zero(double *d, size_t n);

/// create matrix rXc
/// \param r
/// \param c
/// \return
extern double *matrix(size_t r, size_t c);

/// create zero-matrix rXc
/// \param r
/// \param c
/// \return
extern double *zeros(size_t r, size_t c);

/// create identity matrix rXc
/// \param r
/// \return
extern double *identity(size_t r);

/// create vector rX1
/// \param r
/// \return
extern double *vector(size_t r);

/// multiply of quaternion, self-assign unsafe
/// \param p
/// \param q
/// \param m
extern void quaternion_multiply(const double *p, const double *q, double *m);

/// vector dot multiply
/// \param a
/// \param b
/// \param r
/// \return
extern double dot(const double *a, const double *b, size_t r);

/// module of vector
/// \param a
/// \param r
/// \return
extern double norm(double *a, size_t r);

/// normalization of vector
/// \param a
/// \param n
/// \param m
extern void normalize(double *a, size_t n, double *m);

/// cross multiply of vector3, self-assign unsafe
/// \param a
/// \param b
/// \param m
extern void cross3(const double *a, const double *b, double *m);

/// cross multiply matrix of vector3, self-assign unsafe
/// \param a
/// \param m
extern void cross3_matrix(const double *a, double *m);

/// copy of matrix
/// \param a
/// \param r
/// \param c
/// \param m
extern void copy(const double *a, size_t r, size_t c, double *m);

/// matrix plus
/// \param a
/// \param b
/// \param r
/// \param c
/// \param m
extern void plus(const double *a, const double *b, size_t r, size_t c, double *m);

/// matrix subtract
/// \param a
/// \param b
/// \param r
/// \param c
/// \param m
extern void subtract(const double *a, const double *b, size_t r, size_t c, double *m);

/// n*I-A
/// \param a
/// \param n
/// \param r
/// \param m
extern void identity_subtract(const double *a, double n, size_t r, double *m);

/// matrix multiply, self-assign unsafe
/// \param a
/// \param b
/// \param r
/// \param n
/// \param c
/// \param m
extern void multiply(const double *a, const double *b, size_t r, size_t n, size_t c, double *m);

/// scalar multiply
/// \param a
/// \param n
/// \param r
/// \param c
/// \param m
extern void multiply_number(const double *a, double n, size_t r, size_t c, double *m);


#endif //GRIZZ_INS_INERTIAL_MATH_H
