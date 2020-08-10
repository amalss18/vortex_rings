# cython: language_level=3, cdivision=True

import numpy as np
cimport numpy as np
from libc.math cimport pi, sqrt
from cython.parallel import prange

cdef double norm(double a, double b, double c) with gil:
    return sqrt(a*a+b*b+c*c)

cdef double cross_0(double a, double b, double c, double i, double j, double k) with gil:
    cdef double product[3]
    product[0] = b*k - j*c
    return product[0]

cdef double cross_1(double a, double b, double c, double i, double j, double k) with gil:
    cdef double product[3]
    product[1] = c*i - a*k
    return product[1]

cdef double cross_2(double a, double b, double c, double i, double j, double k) with gil:
    cdef double product[3]
    product[2] = a*j - b*i
    return product[2]

cpdef vel_calc(double[:,:] vorticity, double[:] x, double[:] y, double[:] z, double [:,:] velocity):
    cdef int i, j, n, k, l
    n = len(x)
    print(n)
    DEF N = 57760
    i = 0
    cdef double vel_sum[3], ref_particle[3], diff_vec[3]
    cdef double x_new[N], y_new[N], z_new[N], vr_new[N][3], vel_new[N][3]
    cdef double dist, product_0, product_1, product_2
    for k in range(N):
        x_new[k] = x[k]
        y_new[k] = y[k]
        z_new[k] = z[k]
        for l in range(3):
            vr_new[k][l] = vorticity[k][l]
            vel_new[k][l] = velocity[k][l]

    for i in prange(n, nogil=True):
        vel_sum[0] = 0.0
        vel_sum[1] = 0.0
        vel_sum[2] = 0.0
        j = 0
        for j in range(n):
            if i!=j:
                diff_vec[0] = x_new[i] - x_new[j]
                diff_vec[1] = y_new[i] - y_new[j]
                diff_vec[2] = z_new[i] - z_new[j]

                dist = norm(diff_vec[0], diff_vec[1], diff_vec[2])
                product_0 = cross_0(vr_new[j][0], vr_new[j][1], vr_new[j][2], diff_vec[0], diff_vec[1], diff_vec[2])
                product_1 = cross_1(vr_new[j][0], vr_new[j][1], vr_new[j][2], diff_vec[0], diff_vec[1], diff_vec[2])
                product_2 = cross_2(vr_new[j][0], vr_new[j][1], vr_new[j][2], diff_vec[0], diff_vec[1], diff_vec[2])
                vel_sum[0] =vel_sum[0] + product_0 / (dist*dist*dist)
                vel_sum[1] =vel_sum[1] + product_1 / (dist*dist*dist)
                vel_sum[2] =vel_sum[2] + product_2 / (dist*dist*dist)


        vel_new[i][0] = vel_sum[0] / (4*pi)
        vel_new[i][1] = vel_sum[1] / (4*pi)
        vel_new[i][2] = vel_sum[2] / (4*pi)

    return vel_new
