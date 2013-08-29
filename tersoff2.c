/*
 * tersoff2.c: Tersoff 2 potential.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "md.h"
#include "util.h"
#include "tersoff2.h"
#include "tersoff1.h"
#include "tersoff2_params.h"

#define DOUBLE_CMP_EPS DBL_EPSILON      /* Value from float.h. */

/**
 * Force & energy.
 */

/* tersoff2_energy: Computes energy. */
void tersoff2_energy()
{
    int i, j;
    double rij;           /* Distance between atom i and atom j */
    double Uij;           /* Bond energy between atom i and atom j */
    double U = 0;         /* Total energy. */
    double f_cutoff_val;

    for (i = 0; i < natoms - 1; i++) {
        for (j = i + 1; j < natoms; j++) {
            rij = t2_distance(i, j);
            f_cutoff_val = t2_f_cutoff(i, j, rij);

            if (!iszero(f_cutoff_val)) {
                /* If potential is enough large. */
                Uij = f_cutoff_val * (t2_f_repulsive(i, j, rij) + 
                                      t2_b(i, j) * t2_f_attractive(i, j, rij));
                /*
                In Andersen fortran program:
                Uij = f_cutoff_val * (a(i, j, rij) * f_repulsive(i, j, rij) + 
                                      b(i, j, rij) * f_attractive(i, j, rij));
                */
                
                U += Uij;
            }
        }
    }
}

/**
 * Interaction functions 
 */

/* t2_f_cutoff: Smooth cutoff function. */
inline double t2_f_cutoff(int i, int j, double rij)
{
    if (rij < Rij(i, j)) {
        return 1;
    } else if ((rij > Rij(i, j)) && (rij < Sij(i, j))) {
        return 1 / 2 + 1 / 2 * cos(M_PI * (rij - Rij(i, j)) / 
                                   (Sij(i, j) - Rij(i, j)));
    } else { 
        return 0;
    }
}

/* t2_f_repulsive: Repulsive pair potential. */
inline double t2_f_repulsive(int i, int j, double rij)
{
    return Aij(i, j) * exp(-lambda_ij(i, j) * rij);
}

/* t2_f_attractive: Attractive pair potential. */
inline double t2_f_attractive(int i, int j, double rij)
{
    return -Bij(i, j) * exp(-mu_ij(i, j) * rij);
}

/**
 * Parameters
 */

/* t2_b: */
inline double t2_b(int i, int j)
{
    return chi(i, j) * pow(1 + pow(beta(i), n(i)) * 
                           pow(t2_zeta(i, j), n(i)), -1 / (2 * n(i)));
}

/* t2_zeta: */
inline double t2_zeta(int i, int j)
{
    double rik, f_cutoff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        if ((k != i) && (k != j)) {
            rik = t2_distance(i, k);
            f_cutoff_val = t2_f_cutoff(i, k, rik);

            if (!iszero(f_cutoff_val)) {
                sum += f_cutoff_val * omega_ij(i, k) * t2_g(i, j, k);
            }
        }
    }

    return sum;
}

/* t2_g: */
inline double t2_g(int i, int j, int k)
{
    return 1 + pow(c(i), 2) / pow(d(i), 2) -
           pow(c(i), 2) / (pow(d(i), 2) + pow((h(i) - t2_cos_theta(i,j,k)), 2));
}

/* t2_cos_theta: */
inline double t2_cos_theta(int i, int j, int k)
{
    return
        ((x[j] - x[i]) * (x[k] - x[i]) + 
         (y[j] - y[i]) * (y[k] - y[i]) +
         (z[j] - z[i]) * (z[k] - z[i])) / 
        (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2) + pow(z[j] - z[i], 2)) *
         sqrt(pow(x[k] - x[i], 2) + pow(y[k] - y[i], 2) + pow(z[k] - z[i], 2)));
}

/* 
 * t2_distance: Compute distance between atom i and atom j. 
 *              Implements periodic boundary conditions.
 */
inline double t2_distance(int i, int j)
{
    double dx = x[j] - x[i];
    double dy = y[j] - y[i];
    double dz = z[j] - z[i];

    if (dx > Lx / 2) {
        dx = dx - Lx;
    } else if (dx < -Lx / 2) {
        dx = dx + Lx;
    }

    if (dy > Ly / 2) {
        dy = dy - Ly;
    } else if (dy < -Ly / 2) {
        dy = dy + Ly;
    }

    if (dz > Lz / 2) {
        dz = dz - Lz;
    } else if (dz < -Lz / 2) {
        dz = dz + Lz;
    }

    return sqrt(pow(dx, 2) + pow(dy, 2) + pow(dz, 2));
}

/* iszero: Test floating point x is not a zero. */
inline int iszero(double x)
{
    return (x < DOUBLE_CMP_EPS) && (-x > -DOUBLE_CMP_EPS);
}

