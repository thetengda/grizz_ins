//
// Created by tdhuang on 2022/4/12.
//

#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "inertial_math.h"

extern int is_zero(double d) {
    return fabs(d) < 1e-15;
}

extern int all_zero(double *d, size_t n) {
    size_t i;
    for (i = 0; i < n; ++i) {
        if (!is_zero(d[i]))
            return 0;
    }
    return 1;
}

extern double *matrix(size_t r, size_t c) {
    if (r == 0 || c == 0) {
        printf("invalid matrix dim %zux%zu", r, c);
        exit(ERR_MAT_DIM);
        return NULL;
    }
    double *m = malloc(sizeof(double) * r * c);
    return m;
}

extern double *zeros(size_t r, size_t c) {
    double *m = matrix(r, c);
    memset(m, 0, sizeof(double) * r * c);
    return m;
}

extern double *identity(size_t r) {
    int i;
    double *m = matrix(r, r);
    memset(m, 0, sizeof(double) * r * r);
    for (i = 0; i < r; ++i)
        m[MID(r, i, i)] = 1.0;
    return m;
}

extern double *vector(size_t r) {
    return matrix(r, 1);
}

extern void quaternion_multiply(const double *p, const double *q, double *m) {
    m[0] = p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3];
    m[1] = p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2];
    m[2] = p[0] * q[2] + p[2] * q[0] + p[3] * q[1] - p[1] * q[3];
    m[3] = p[0] * q[3] + p[3] * q[0] + p[1] * q[2] - p[2] * q[1];
}

extern double dot(const double *a, const double *b, size_t r) {
    int i;
    double ret = 0;
    for (i = 0; i < r; ++i)
        ret += a[i] * b[i];
    return ret;
}

extern double norm(double *a, size_t r) {
    return sqrt(dot(a, a, r));
}

extern void normalize(double *a, size_t n, double *m) {
    size_t i;
    double sum = 0;
    if (all_zero(a, n)) {
        printf("vector all 0");
        exit(ERR_NORM_ZERO);
    }

    for (i = 0; i < n; ++i) {
        sum += SQR(a[i]);
    }
    sum = sqrt(sum);
    for (i = 0; i < n; ++i) {
        m[i] = a[i] / sum;
    }
}

extern void cross3(const double *a, const double *b, double *m) {
    m[0] = a[1] * b[2] - a[2] * b[1];
    m[1] = a[2] * b[0] - a[0] * b[2];
    m[2] = a[0] * b[1] - a[1] * b[0];
}

extern void cross3_matrix(const double *a, double *m) {
    m[MID(3, 0, 0)] = 0;
    m[MID(3, 0, 1)] = -a[2];
    m[MID(3, 0, 2)] = a[1];
    m[MID(3, 1, 0)] = a[2];
    m[MID(3, 1, 1)] = 0;
    m[MID(3, 1, 2)] = -a[0];
    m[MID(3, 2, 0)] = -a[1];
    m[MID(3, 2, 1)] = a[0];
    m[MID(3, 2, 2)] = 0;
}

extern void copy(const double *a, size_t r, size_t c, double *m) {
    memmove(m, a, sizeof(double) * r * c);
}

extern void plus(const double *a, const double *b, size_t r, size_t c, double *m) {
    size_t i, j;

    for (i = 0; i < r; ++i)
        for (j = 0; j < c; ++j)
            m[MID(c, i, j)] = a[MID(c, i, j)] + b[MID(c, i, j)];

}

extern void subtract(const double *a, const double *b, size_t r, size_t c, double *m) {
    size_t i, j;

    for (i = 0; i < r; ++i)
        for (j = 0; j < c; ++j)
            m[MID(c, i, j)] = a[MID(c, i, j)] - b[MID(c, i, j)];

}

extern void identity_subtract(const double *a, double n, size_t r, double *m) {
    size_t i, j;
    for (i = 0; i < r; ++i) {
        for (j = 0; j < r; ++j) {
            m[MID(r, i, j)] = (i == j ? n : 0) - a[MID(r, i, j)];
        }
    }

}

extern void multiply(const double *a, const double *b, size_t r, size_t n, size_t c, double *m) {
    size_t i, j, k;
    for (i = 0; i < r; ++i) {
        for (j = 0; j < c; ++j) {
            m[MID(c, i, j)] = 0;
            for (k = 0; k < n; ++k) {
                m[MID(c, i, j)] += a[MID(n, i, k)] * b[MID(c, k, j)];
            }
        }
    }
}

extern void multiply_number(const double *a, double n, size_t r, size_t c, double *m) {
    size_t i, j;

    for (i = 0; i < r; ++i) {
        for (j = 0; j < c; ++j) {
            m[MID(c, i, j)] = a[MID(c, i, j)] * n;
        }
    }

}


