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

static double Uij_diff         (atom_t *az, atom_t *ai, atom_t *aj, 
                                coord_t coord);
static double f_cutoff_diff    (atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord); 
static double rij_diff         (atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord);
static double f_attractive_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord);
static double f_repulsive_diff (atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord);
static double bij_diff         (atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord);
static double zeta_diff        (atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord);
static double g_diff           (atom_t *az, atom_t *ai, atom_t *aj, atom_t *ak, 
                                double rij, double rik, coord_t coord);
static double cos_theta_diff   (atom_t *az, atom_t *ai, atom_t *aj, atom_t *ak, 
                                double rij, double rik, coord_t coord);
static int kronecker_delta(atom_t *ai, atom_t *aj);
static double getcoord(atom_t *ai, coord_t coord);

/* tersoff2_forces: Computes forces. */
void tersoff2_forces()
{
    int z, i, j;

    for (z = 0; z < natoms; z++) {
        double Fx = 0, Fy = 0, Fz = 0;  /* Forces */
        atom_t *az = &atom[z];
        atom_t *ai = &atom[z];
        atom_t *aj = &atom[z];

        printf("z = %d\n", z);
        Fx = Fy = Fz = 0;

        /* Compute sum of differential of U. */
        for (i = 0; i < natoms; i++) {
            printf("i = %d\n", i);

            for (j = 0; j < natoms; j++) {
                if (i != j) { 
                    Fx += Uij_diff(az, ai, aj, COORD_X);
                    Fy += Uij_diff(az, ai, aj, COORD_Y);
                    Fz += Uij_diff(az, ai, aj, COORD_Z);
                }
            }

            /*
            printf("Fx = %f\n", Fx);
            printf("Fy = %f\n", Fy);
            printf("Fz = %f\n", Fz);
            */
        }

        Fx *= -0.5;
        Fy *= -0.5;
        Fz *= -0.5;

        atom[z].ax = Fx / mass(az);
        atom[z].ay = Fy / mass(az);
        atom[z].az = Fz / mass(az);
    }
}

/* Uij_diff: */
static double Uij_diff(atom_t *az, atom_t *ai, atom_t *aj, coord_t coord)
{
    double rij, f_cutoff_val, f_cutoff_diff_val, Uij_diff_val = 0; 

    rij = t2_distance(ai, aj);

    /* If f_cutoff = 0 then f_cutoff_diff = 0 */

    f_cutoff_val = t2_f_cutoff(ai, aj, rij);

    if (!iszero(f_cutoff_val)) {

        Uij_diff_val = f_cutoff_val * 
                       (f_repulsive_diff(az, ai, aj, rij, coord) + 
                        bij_diff(az, ai, aj, rij, coord) * 
                        t2_f_attractive(ai, aj, rij) +
                        t2_b(ai, aj) * 
                        f_attractive_diff(az, ai, aj, rij, coord));

        f_cutoff_diff_val = f_cutoff_diff(az, ai, aj, rij, coord);

        if (!iszero(f_cutoff_diff_val)) {
            Uij_diff_val += f_cutoff_diff_val * 
                            (t2_f_repulsive(ai, aj, rij) + 
                             t2_b(ai, aj) * t2_f_attractive(ai, aj, rij));
        }

    }

    return Uij_diff_val;
}

/* f_cutoff_diff: */
static double f_cutoff_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                            coord_t coord)
{
    if (rij < Rij(ai, aj)) {
        return 0;
    } else if ((rij > Rij(ai, aj)) && (rij < Sij(ai, aj))) {
        return - M_PI / (2 * (Sij(ai, aj), Rij(ai, aj))) *
                 cos(M_PI * (rij - Rij(ai, aj)) / (Sij(ai, aj) - Rij(ai, aj))) *
                 rij_diff(az, ai, aj, rij, coord);
    } else /* if (rij > Sij(i, j)) */ {
        return 0;
    }
}

/* rij_diff: */
static double rij_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                       coord_t coord)
{
    if (az == ai) {
        return (getcoord(ai, coord) - getcoord(aj, coord)) / rij;
    } else if (az == aj) {
        return (getcoord(aj, coord) - getcoord(ai, coord)) / rij;
    } else {
        return 0;
    }
}

/* f_attractive_diff: */
static double f_attractive_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                                coord_t coord)
{
    return -lambda_ij(ai, aj) * Aij(ai, aj) * exp(-lambda_ij(ai, aj) * rij) *
            rij_diff(az, ai, aj, rij, coord);
}

/* f_repulsive_diff: */
static double f_repulsive_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                               coord_t coord)
{
    return -lambda_ij(ai, aj) * -Bij(ai, aj) * exp(-lambda_ij(ai, aj) * rij) *
            rij_diff(az, ai, aj, rij, coord);
}

/* bij_diff: */
static double bij_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                       coord_t coord)
{
    double bij_denom_tmp, t2_zeta_val;

    /* 
     * Fork in bij_diff implementation see Andersen:
     * force_tersoff.f90:193-200
     */

    bij_denom_tmp = 1 + pow(beta(ai), n(ai)) * pow(t2_zeta(ai, aj), n(ai));

    if (iszero(bij_denom_tmp)) {
        return 0.0;
    } else {
        t2_zeta_val = t2_zeta(ai, aj);

        if (iszero(t2_zeta_val)) {
            return 0.0;
        } else {
            return chi(ai, aj) * (-pow(beta(ai), n(ai)) * 
                                 pow(t2_zeta_val, n(ai) - 1)) / 
                   (2 * pow(bij_denom_tmp, 1 / (2 * n(ai)) + 1)) * 
                   zeta_diff(az, ai, aj, rij, coord);
        }
    }
}

/* zeta_diff: */
static double zeta_diff(atom_t *az, atom_t *ai, atom_t *aj, double rij, 
                        coord_t coord)
{
    double rik, f_cutoff_val, f_cutoff_diff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        atom_t *ak = &atom[k];

        if ((ak != ai) && (ak != aj)) {
            rik = t2_distance(ai, ak);

            /* If f_cutoff = 0 then f_cutoff_diff = 0 */

            f_cutoff_val = t2_f_cutoff(ai, ak, rik);

            if (!iszero(f_cutoff_val)) {

                sum = f_cutoff_val * omega_ij(ai, ak) * 
                      g_diff(az, ai, aj, ak, rij, rik, coord);

                f_cutoff_diff_val = f_cutoff_diff(az, ai, ak, rik, coord);

                if (!iszero(f_cutoff_diff_val)) {
                    sum += f_cutoff_diff_val * omega_ij(ai, ak) * 
                           t2_g(ai, aj, ak);
                }
            }
        }
    }

    return sum;
}

/* g_diff: */
static double g_diff(atom_t *az, atom_t *ai, atom_t *aj, atom_t *ak, 
                     double rij, double rik, coord_t coord)
{
    return -(2 * pow(c(ai), 2) * (h(ai) - t2_cos_theta(ai, aj, ak))) / 
            pow(pow(d(ai), 2) + 
                pow(h(ai) - t2_cos_theta(ai, aj, ak), 2), 2) * 
            cos_theta_diff(az, ai, aj, ak, rij, rik, coord);
}

/* cos_theta_diff: */
static double cos_theta_diff(atom_t *az, atom_t *ai, atom_t *aj, atom_t *ak, 
                             double rij, double rik, coord_t coord)
{
    double i_coord = getcoord(ai, coord);
    double j_coord = getcoord(aj, coord);
    double k_coord = getcoord(ak, coord);

    double rij_vec = j_coord - i_coord;
    double rik_vec = k_coord - i_coord;

    int delta_zk = kronecker_delta(az, ak);
    int delta_zi = kronecker_delta(az, ai);
    int delta_zj = kronecker_delta(az, aj);

    return (rij_vec * (delta_zk - delta_zi) + 
            rik_vec * (delta_zj - delta_zi)) / (rij * rik) - 
           t2_cos_theta(ai, aj, ak) * 
           ((rij_vec / pow(rij, 2)) * (delta_zj - delta_zi) +
            (rik_vec / pow(rik, 2)) * (delta_zk - delta_zi));
}

/* kronecker_delta: */ 
static int kronecker_delta(atom_t *ai, atom_t *aj)
{
    return ai == aj ? 1 : 0;
}

/* 
 * getcoord: Return coordinate x, y or z of atom i 
 *            according to the coord variable.
 */
static double getcoord(atom_t *ai, coord_t coord)
{
    switch(coord) {
        case COORD_X:
            return ai->x;
            break;
        case COORD_Y:
            return ai->y;
            break;
        case COORD_Z:
            return ai->z;
            break;
    }
    return 0;
}
