/*
 * tersoff2.c: Tersoff 2 potential - forces.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "md.h"
#include "util.h"
#include "tersoff2.h"
#include "tersoff2_forces.h"
#include "tersoff2_params.h"
#include "celllist.h"

typedef enum {
    COORD_X,
    COORD_Y,
    COORD_Z
} coord_t;

/*
 * z     - is the current atom,
 * i, j  - neighbour atoms (to z),
 * rij   - distance between atoms i and j,
 * coord - x, y or z - coordinates
 */

/* This variables are global to pass into celllist_scan_atoms. */
double Fx = 0, Fy = 0, Fz = 0;  /* Forces */
int atom_z = 0;                 /* Current atom */

static void tersoff2_forces_ij(int atom_i, int atom_j);

static double Uij_diff         (int z, int i, int j, coord_t coord);
static double f_cutoff_diff    (int z, int i, int j, double rij, coord_t coord);
static double rij_diff         (int z, int i, int j, double rij, coord_t coord);
static double f_attractive_diff(int z, int i, int j, double rij, coord_t coord);
static double f_repulsive_diff (int z, int i, int j, double rij, coord_t coord);
static double bij_diff         (int z, int i, int j, double rij, coord_t coord);
static double zeta_diff        (int z, int i, int j, double rij, coord_t coord);
static double g_diff           (int z, int i, int j, int k, 
                                double rij, double rik, coord_t coord);
static double cos_theta_diff   (int z, int i, int j, int k, 
                                double rij, double rik, coord_t coord);
static int kronecker_delta(int i, int j);
static double getcoord(int i, coord_t coord);

/* tersoff2_forces: Computes forces. */
void tersoff2_forces()
{

    /* Main loop by atoms: compute force for each atom. */
    for (atom_z = 0; atom_z < natoms; atom_z++) {

        printf("z = %d\n", atom_z);
        Fx = Fy = Fz = 0;

        celllist_scan_atoms(tersoff2_forces_ij);

        Fx *= -0.5;
        Fy *= -0.5;
        Fz *= -0.5;

        ax[atom_z] = Fx / mass(atom_z);
        ay[atom_z] = Fy / mass(atom_z);
        az[atom_z] = Fz / mass(atom_z);
    }
}

/* tersoff2_energy_ij: Compute energy for atom pair (atom_i, atom_j). */
static void tersoff2_forces_ij(int atom_i, int atom_j)
{
    if (atom_i != atom_j) { 
        Fx += Uij_diff(atom_z, atom_i, atom_j, COORD_X);
        Fy += Uij_diff(atom_z, atom_i, atom_j, COORD_Y);
        Fz += Uij_diff(atom_z, atom_i, atom_j, COORD_Z);

        /*     printf("Fx = %f\n", Fx); */
        /*     printf("Fy = %f\n", Fy); */
        /*     printf("Fz = %f\n", Fz); */
    }
}

/* Uij_diff: */
static double Uij_diff(int z, int i, int j, coord_t coord)
{
    double rij, f_cutoff_val, f_cutoff_diff_val, Uij_diff_val = 0; 

    rij = t2_distance(i, j);

    /* If f_cutoff = 0 then f_cutoff_diff = 0 */

    f_cutoff_val = t2_f_cutoff(i, j, rij);

    if (!iszero(f_cutoff_val)) {

        Uij_diff_val = f_cutoff_val * 
                       (f_repulsive_diff(z, i, j, rij, coord) + 
                        bij_diff(z, i, j, rij, coord) * 
                        t2_f_attractive(i, j, rij) +
                        t2_b(i, j) * 
                        f_attractive_diff(z, i, j, rij, coord));

        f_cutoff_diff_val = f_cutoff_diff(z, i, j, rij, coord);

        if (!iszero(f_cutoff_diff_val)) {
            Uij_diff_val += f_cutoff_diff_val * 
                            (t2_f_repulsive(i, j, rij) + 
                             t2_b(i, j) * t2_f_attractive(i, j, rij));
        }

    }

    return Uij_diff_val;
}

/* f_cutoff_diff: */
static double f_cutoff_diff(int z, int i, int j, double rij, coord_t coord)
{
    if (rij < Rij(i, j)) {
        return 0;
    } else if ((rij > Rij(i, j)) && (rij < Sij(i, j))) {
        return - M_PI / (2 * (Sij(i, j), Rij(i, j))) *
                 cos(M_PI * (rij - Rij(i, j)) / (Sij(i, j) - Rij(i, j))) *
                 rij_diff(z, i, j, rij, coord);
    } else /* if (rij > Sij(i, j)) */ {
        return 0;
    }
}

/* rij_diff: */
static double rij_diff(int z, int i, int j, double rij, coord_t coord)
{
    if (z == i) {
        return (getcoord(i, coord) - getcoord(j, coord)) / rij;
    } else if (z == j) {
        return (getcoord(j, coord) - getcoord(i, coord)) / rij;
    } else {
        return 0;
    }
}

/* f_attractive_diff: */
static double f_attractive_diff(int z, int i, int j, double rij, coord_t coord)
{
    return -lambda_ij(i, j) * Aij(i, j) * exp(-lambda_ij(i, j) * rij) *
            rij_diff(z, i, j, rij, coord);
}

/* f_repulsive_diff: */
static double f_repulsive_diff(int z, int i, int j, double rij, coord_t coord)
{
    return -lambda_ij(i, j) * -Bij(i, j) * exp(-lambda_ij(i, j) * rij) *
            rij_diff(z, i, j, rij, coord);
}

/* bij_diff: */
static double bij_diff(int z, int i, int j, double rij, coord_t coord)
{
    double bij_denom_tmp, t2_zeta_val;

    /* 
     * Fork in bij_diff implementation see Andersen:
     * force_tersoff.f90:193-200
     */

    bij_denom_tmp = 1 + pow(beta(i), n(i)) * pow(t2_zeta(i, j), n(i));

    if (iszero(bij_denom_tmp)) {
        return 0.0;
    } else {
        t2_zeta_val = t2_zeta(i, j);

        if (iszero(t2_zeta_val)) {
            return 0.0;
        } else {
            return chi(i, j) * (-pow(beta(i), n(i)) * 
                                 pow(t2_zeta_val, n(i) - 1)) / 
                   (2 * pow(bij_denom_tmp, 1 / (2 * n(i)) + 1)) * 
                   zeta_diff(z, i, j, rij, coord);
        }
    }
}

/* zeta_diff: */
static double zeta_diff(int z, int i, int j, double rij, coord_t coord)
{
    double rik, f_cutoff_val, f_cutoff_diff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        if ((k != i) && (k != j)) {
            rik = t2_distance(i, k);

            /* If f_cutoff = 0 then f_cutoff_diff = 0 */

            f_cutoff_val = t2_f_cutoff(i, k, rik);

            if (!iszero(f_cutoff_val)) {

                sum = f_cutoff_val * omega_ij(i, k) * 
                      g_diff(z, i, j, k, rij, rik, coord);

                f_cutoff_diff_val = f_cutoff_diff(z, i, k, rik, coord);

                if (!iszero(f_cutoff_diff_val)) {
                    sum += f_cutoff_diff_val * omega_ij(i, k) * t2_g(i, j, k);
                }
            }
        }
    }

    return sum;
}

/* g_diff: */
static double g_diff(int z, int i, int j, int k, 
                     double rij, double rik, coord_t coord)
{
    return -(2 * pow(c(i), 2) * (h(i) - t2_cos_theta(i, j, k))) / 
            pow(pow(d(i), 2) + 
                pow(h(i) - t2_cos_theta(i, j, k), 2), 2) * 
            cos_theta_diff(z, i, j, k, rij, rik, coord);
}

/* cos_theta_diff: */
static double cos_theta_diff(int z, int i, int j, int k, 
                             double rij, double rik, coord_t coord)
{
    double i_coord = getcoord(i, coord);
    double j_coord = getcoord(j, coord);
    double k_coord = getcoord(k, coord);

    double rij_vec = j_coord - i_coord;
    double rik_vec = k_coord - i_coord;

    int delta_zk = kronecker_delta(z, k);
    int delta_zi = kronecker_delta(z, i);
    int delta_zj = kronecker_delta(z, j);

    return (rij_vec * (delta_zk - delta_zi) + 
            rik_vec * (delta_zj - delta_zi)) / (rij * rik) - 
           t2_cos_theta(i, j, k) * 
           ((rij_vec / pow(rij, 2)) * (delta_zj - delta_zi) +
            (rik_vec / pow(rik, 2)) * (delta_zk - delta_zi));
}

/* kronecker_delta: */ 
static int kronecker_delta(int i, int j)
{
    return i == j ? 1 : 0;
}

/* 
 * getcoord: Return coordinate x, y or z of atom i 
 *            according to the coord variable.
 */
static double getcoord(int i, coord_t coord)
{
    switch(coord) {
        case COORD_X:
            return x[i];
            break;
        case COORD_Y:
            return y[i];
            break;
        case COORD_Z:
            return z[i];
            break;
    }
    return 0;
}
