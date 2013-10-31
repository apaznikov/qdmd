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
#include "celllist.h"

#define DOUBLE_CMP_EPS DBL_EPSILON      /* Value from float.h. */

/**
 * Force & energy.
 */

double U = 0;         /* Total energy. */

static void tersoff2_energy_ij(int atom_i, int atom_j);

/* tersoff2_energy: Computes energy. */
void tersoff2_energy()
{
    celllist_scan_atoms(tersoff2_energy_ij);

    /* FIXME: for DEBUG */
    exit(0);
}

/* tersoff2_energy_ij: Compute energy for atom pair (atom_i, atom_j). */
static void tersoff2_energy_ij(int atom_i, int atom_j)
{
    double rij;           /* Distance between atom i and atom j */
    double f_cutoff_val;

    /* Compute energy */
    rij = t2_distance_noPBC(atom_i, atom_j);
    f_cutoff_val = t2_f_cutoff(atom_i, atom_j, rij);

    if (!iszero(f_cutoff_val)) {
        /* Bond energy between atom i and atom j */
        double Uij;           
        /* If potential is enough large. */
        Uij = f_cutoff_val * 
              (t2_f_repulsive(atom_i, atom_j, rij) + 
               t2_b(atom_i, atom_j) * 
               t2_f_attractive(atom_i, atom_j, rij));

        /*
        In Andersen fortran program:
        Uij = f_cutoff_val * 
              (a(atom_i, atom_j, rij) * 
               f_repulsive(atom_i, atom_j, rij) + 
               b(atom_i, atom_j, rij) * 
               f_attractive(atom_i, atom_j, rij));
        */
        
        U += Uij;
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
inline double t2_distance(int atom_i, int atom_j)
{
    double dx = x[atom_j] - x[atom_i];
    double dy = y[atom_j] - y[atom_i];
    double dz = z[atom_j] - z[atom_i];

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

/* 
 * t2_distance_noPBC: Compute distance between atom i and atom j. 
 *                    NOT implements periodic boundary conditions.
 */
inline double t2_distance_noPBC(int atom_i, int atom_j)
{
    return sqrt(pow(x[atom_j] - x[atom_i], 2) + 
                pow(y[atom_j] - y[atom_i], 2) + 
                pow(z[atom_j] - z[atom_i], 2));
}

/* iszero: Test floating point x is not a zero. */
inline int iszero(double x)
{
    return (x < DOUBLE_CMP_EPS) && (-x > -DOUBLE_CMP_EPS);
}

