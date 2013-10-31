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

/* tersoff2_energy: Computes energy. */
void tersoff2_energy()
{
    double U = 0;         /* Total energy. */
    int n_int = 0;

    cell_t cell;
    cell_t neigh_cell_raw;  /* Raw cell for scan all neighbour cells. */
    cell_t neigh_cell;      /* Working neighbour cell. */

    init_cell(&cell);

    while (scan_cells(&cell)) {
        cell_vec_to_scal(&cell);

        init_neigh_cell(cell, &neigh_cell_raw);
        while (scan_neigh_cells(cell, &neigh_cell_raw)) {
            int atom_cell, atom_neighcell;
            double atom_shift_x, atom_shift_y, atom_shift_z;

            cell_periodic_bound_cond(neigh_cell_raw, &neigh_cell,
                                     &atom_shift_x, &atom_shift_y, 
                                     &atom_shift_z);
            cell_vec_to_scal(&neigh_cell);

            /* Now scan over the atoms in cell and atoms in neighbour cell. */

            if (!init_atom_in_cell(cell, &atom_cell)) {
                continue;   /* if no atoms in the cell */
            }

            do {
                if (!init_atom_in_cell(neigh_cell, &atom_neighcell)) {
                    continue;   /* if no atoms in the cell */
                }

                do {

                    /* Avoid double counting of pair (i, j) */
                    if (atom_cell < atom_neighcell) {

                        /* TODO implement PBC via shift vector. */
                        atom_t *ai = &atom[atom_cell];
                        double rij, f_cutoff_val;

                        /* 
                         * Atomic positions in an outside cell must be treated
                         * as they are and must not be pulled back into the
                         * central simuation box.
                         */
                        atom_t aj = atom[atom_neighcell];
                        aj.x += atom_shift_x;
                        aj.y += atom_shift_y;
                        aj.z += atom_shift_z;

                        /* Compute energy */
                        rij = t2_distance_noPBC(ai, &aj); 
                        /* rij = t2_distance(ai, &aj); */
                        f_cutoff_val = t2_f_cutoff(ai, &aj, rij);

                        if (!iszero(f_cutoff_val)) {
                            /* Bond energy between atom i and atom j */
                            double Uij;           

                            n_int++;
                            /* printf("n_int = %d\n", n_int); */

                            /* If potential is enough large. */
                            Uij = f_cutoff_val * 
                                  (t2_f_repulsive(ai, &aj, rij) + 
                                   t2_b(ai, &aj) * t2_f_attractive(ai, &aj, rij));

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
                } while (scan_atom_in_cell(neigh_cell, &atom_neighcell));
            } while (scan_atom_in_cell(cell, &atom_cell));
        }
    }

    printf("Total energy U = %f\n", U);
    printf("n_int = %d\n", n_int);

    exit(0);
}

/**
 * Interaction functions 
 */

/* t2_f_cutoff: Smooth cutoff function. */
inline double t2_f_cutoff(atom_t *ai, atom_t *aj, double rij)
{
    if (rij < Rij(ai, aj)) {
        return 1;
    } else if ((rij > Rij(ai, aj)) && (rij < Sij(ai, aj))) {
        return 1 / 2 + 1 / 2 * cos(M_PI * (rij - Rij(ai, aj)) / 
                                   (Sij(ai, aj) - Rij(ai, aj)));
    } else { 
        return 0;
    }
}

/* t2_f_repulsive: Repulsive pair potential. */
inline double t2_f_repulsive(atom_t *ai, atom_t *aj, double rij)
{
    return Aij(ai, aj) * exp(-lambda_ij(ai, aj) * rij);
}

/* t2_f_attractive: Attractive pair potential. */
inline double t2_f_attractive(atom_t *ai, atom_t *aj, double rij)
{
    return -Bij(ai, aj) * exp(-mu_ij(ai, aj) * rij);
}

/**
 * Parameters
 */

/* t2_b: */
inline double t2_b(atom_t *ai, atom_t *aj)
{
    return chi(ai, aj) * pow(1 + pow(beta(ai), n(ai)) * 
                           pow(t2_zeta(ai, aj), n(ai)), -1 / (2 * n(ai)));
}

/* t2_zeta: */
inline double t2_zeta(atom_t *ai, atom_t *aj)
{
    double rik, f_cutoff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        atom_t *ak = &atom[k];

        if ((ak != ai) && (ak != aj)) {
            rik = t2_distance(ai, ak);
            f_cutoff_val = t2_f_cutoff(ai, ak, rik);

            if (!iszero(f_cutoff_val)) {
                sum += f_cutoff_val * omega_ij(ai, ak) * t2_g(ai, aj, ak);
            }
        }
    }

    return sum;
}

/* t2_g: */
inline double t2_g(atom_t *ai, atom_t *aj, atom_t *ak)
{
    return 1 + pow(c(ai), 2) / pow(d(ai), 2) -
           pow(c(ai), 2) / (pow(d(ai), 2) + pow((h(ai) - 
                            t2_cos_theta(ai, aj, ak)), 2));
}

/* t2_cos_theta: */
inline double t2_cos_theta(atom_t *ai, atom_t *aj, atom_t *ak)
{
    return
        ((aj->x - ai->x) * (ak->x - ai->x) + 
         (aj->y - ai->y) * (ak->y - ai->y) +
         (aj->z - ai->z) * (ak->z - ai->z)) / 
        (sqrt(pow(aj->x - ai->x, 2) + 
              pow(aj->y - ai->y, 2) + 
              pow(aj->z - ai->z, 2)) *
         sqrt(pow(ak->x - ai->x, 2) + 
              pow(ak->y - ai->y, 2) + 
              pow(ak->z - ai->z, 2)));
}

/* 
 * t2_distance: Compute distance between atom i and atom j. 
 *              Implements periodic boundary conditions.
 */
inline double t2_distance(atom_t *ai, atom_t *aj)
{
    double dx = aj->x - ai->x;
    double dy = aj->y - ai->y;
    double dz = aj->z - ai->z;

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
inline double t2_distance_noPBC(atom_t *ai, atom_t *aj)
{
    return sqrt(pow(aj->x - ai->x, 2) + 
                pow(aj->y - ai->y, 2) + 
                pow(aj->z - ai->z, 2));
}

/* iszero: Test floating point x is not a zero. */
inline int iszero(double x)
{
    return (x < DOUBLE_CMP_EPS) && (-x > -DOUBLE_CMP_EPS);
}
