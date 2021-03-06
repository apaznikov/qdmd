/*
 * cellist.c: Linked cell list implementation
 * 
 * Copyright (C) 2013 Alexey Paznikov
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#include "md.h"
#include "util.h"
#include "celllist.h"

/* End of a list flag. */
enum { LIST_END = -1 }; 

/* Addition to r_cutoff, this value from [Rapaport, 55, ch. 3.4] */
static const double delta_r = 0.35; 

/* 
 * Cutoff distance: 
 * Larger than this disance interactions can be neglected.
 * 
 * r_cutoff = 2.5: [Rapaport, 48], [Allen, Tildesley, 148] 
 */
static const double r_cutoff = 2.5;        

/* Number of cells by each coordinate. */
/* NOTE: possibly double not int? */
static int cell_Lx;
static int cell_Ly;
static int cell_Lz;
static int ncells;
static int cell_LyLz;    /* temp */

/* Size of a one cell. */
static double cell_rx;
static double cell_ry;
static double cell_rz;

/* Linked list: head and list as such. */
static int *list, *head;

static void init_cell(cell_t *cell);
static bool scan_cells(cell_t *cell);
static void init_neigh_cell(cell_t cell,  cell_t *neigh_cell);
static bool scan_neigh_cells(cell_t cell, cell_t *neigh_cell);

static bool init_atom_in_cell(cell_t cell, int *atom_i);
static bool scan_atom_in_cell(cell_t cell, int *atom_i);

static void cell_vec_to_scal(cell_t *cell);
static void cell_periodic_bound_cond(cell_t cell_orig, cell_t *cell_new,
                              double *atom_shift_x, double *atom_shift_y, 
                              double *atom_shift_z);

/* 
 * celllist_scan_atoms: Scan over all atom pairs. 
 *                      f - function (forces or energy) to do with pair of atoms
 */
void celllist_scan_atoms(void (*func) (int atom_i, int atom_j))
{
    cell_t cell;
    cell_t neigh_cell_raw;  /* Cell for scan all neighbour cells without PBC. */
    cell_t neigh_cell;      /* Working neighbour cell. */

    init_cell(&cell);

    while (scan_cells(&cell)) {

        cell_vec_to_scal(&cell);
        if (cell.scal % 1000 == 0) {
            printf("cell = %d\n", cell.scal);
        }

        init_neigh_cell(cell, &neigh_cell_raw);

        while (scan_neigh_cells(cell, &neigh_cell_raw)) {

            int atom_cell, atom_neighcell;
            double atom_shift_x, atom_shift_y, atom_shift_z;

            cell_periodic_bound_cond(neigh_cell_raw, &neigh_cell,
                                   &atom_shift_x, &atom_shift_y, &atom_shift_z);
            cell_vec_to_scal(&neigh_cell);

            if (cell.scal >= neigh_cell.scal) {
                continue;   /* avoid double counting of cell pairs */
            }

            /* Now scan over the atoms in cell and atoms in neighbour cell. */

            if (!init_atom_in_cell(cell, &atom_cell)) {
                continue;   /* if no atoms in the cell */
            }

            do {
                if (!init_atom_in_cell(neigh_cell, &atom_neighcell)) {
                    continue;   /* if no atoms in the cell */
                }

                do {
                    /* fprintf(stderr, "cell %d neigh %d a_i %d a_j %d\n", */
                    /* cell.scal, neigh_cell.scal, atom_i, atom_j); */

                    int atom_i = atom_cell;
                    int atom_j = atom_neighcell;

                    /* Avoid double counting of pair (i, j) */
                    if (atom_i >= atom_j) {
                        continue;
                    }

                    /* 
                     * TODO 
                     * implement PBC via shift vector. 
                     */

                    /* 
                     * TODO
                     * Atomic positions in an outside cell must be treated
                     * as they are and must not be pulled back into the 
                     * central simuation box.
                     */
                    /*
                    atom_t aj = atom[atom_neighcell];
                    aj.x += atom_shift_x;
                    aj.y += atom_shift_y;
                    aj.z += atom_shift_z;
                    */
                    
                    /* 
                     * Call function to compute energy (force, etc) 
                     * for atom pair.
                     */
                    func(atom_i, atom_j);

                } while (scan_atom_in_cell(neigh_cell, &atom_neighcell));
            } while (scan_atom_in_cell(cell, &atom_cell));
        }
    }
}

/* linked_cell_init: Init head and all constants. */
void linked_cell_init()
{
    /* 
     * r_cutoff_more is more than simply r_cutoff in order to don't 
     * update cells after each step. 
     */
    const double r_cutoff_more = r_cutoff + delta_r;

    int cell;

    /* Size of lattice: number of cells by coordinates */
    cell_Lx = floor(Lx / r_cutoff_more);
    cell_Ly = floor(Ly / r_cutoff_more);
    cell_Lz = floor(Lz / r_cutoff_more);

    ncells = cell_Lx * cell_Ly * cell_Lz;

    /* temp */
    cell_LyLz = cell_Ly * cell_Lz;   

    /* Size of a one cell - more than r_cutoff_more. */
    cell_rx = Lx / cell_Lx;
    cell_ry = Ly / cell_Ly;
    cell_rz = Lz / cell_Lz;

    list = xmalloc(sizeof(*list) * natoms);
    head = xmalloc(sizeof(*head) * ncells);

    /* Reset headers to emtpy (end) */
    for (cell = 0; cell < ncells; cell++) {
        head[cell] = LIST_END;
    }
}

/* linked_cell_update: Build list of atoms belongs to each cell. */
void linked_cell_update()
{
    int atom;
    cell_t cell;

    /* Construct head and list */
    for (atom = 0; atom < natoms; atom++) {

        /* Cell index (vector) to which this atom belongs */
        cell.x = x[atom] / cell_rx;
        cell.y = y[atom] / cell_ry;
        cell.z = z[atom] / cell_rz;

        cell_vec_to_scal(&cell);

        /* Link to the previous occupant */
        list[atom] = head[cell.scal];

        /* The last one goes to the header */
        head[cell.scal] = atom;
    }
}

/* cell_vec_to_scal: Convert coordinates of cell from vector to scalar. */
static void cell_vec_to_scal(cell_t *cell)
{
    cell->scal = cell->x * cell_LyLz + cell->y * cell_Lz + cell->z;
}

/* linked_cell_finalize(): */
void linked_cell_finalize()
{
    free(list);
    free(head);
}

/* 
 * cell_update_test: Test whether the time for cell list update has come. 
 *                   This criterion was taken from [Rapoport, 55, chapter 3.4] 
 */
bool cell_update_test()
{
    int atom;
    double max_v = vx[0];
    static double sum_max_v = 0;

    /* Find maximum of velocities. */
    for (atom = 0; atom < natoms; atom++) {

        if (vx[atom] > max_v) {
            max_v = vx[atom];
        }

        if (vy[atom] > max_v) {
            max_v = vy[atom];
        }

        if (vz[atom] > max_v) {
            max_v = vz[atom];
        }
    }

    /* Compute sum of velocities. */
    sum_max_v += max_v;

    /* This criterion was taken from [Rapoport, 55, chapter 3.4] */
    if (sum_max_v > delta_r / (2 * timestep_dt)) {
        return true;
        sum_max_v = 0;
    } else {
        return false;
    }

    return true;
}

/* init_cell: Before start scan, initialize cell. */
static void init_cell(cell_t *cell)
{
    cell->x = 0;
    cell->y = 0;
    cell->z = -1; /* -1 for first increment in scan_cells */
}

/* scan_cells: Scan over cells. Goes to the next cell. */
static bool scan_cells(cell_t *cell)
{
    /* Increment the cell coordinates. Check bounds. */
    cell->z++;
    if (cell->z == cell_Lz) {
        cell->z = 0;
        cell->y++;

        if (cell->y == cell_Ly) {
            cell->y = 0;
            cell->x++;

            if (cell->x == cell_Lx) {
                cell->x = 0;
                return false;
            }
        }
    }

    return true;
}

/* init_neigh_cell: Before start scan, initialize cell. */
static void init_neigh_cell(cell_t cell, cell_t *neigh_cell)
{
    neigh_cell->x = cell.x - 1;
    neigh_cell->y = cell.y - 1;
    neigh_cell->z = cell.z - 1 - 1; /* -1 for first increment */
}

/* scan_neigh_cells: */
static bool scan_neigh_cells(cell_t cell, cell_t *neigh_cell)
{
    /* Increment cell coordinates. Check bounds. */
    neigh_cell->z++;
    if (neigh_cell->z > cell.z + 1) {
        neigh_cell->z = cell.z - 1;
        neigh_cell->y++;

        if (neigh_cell->y > cell.y + 1) {
            neigh_cell->y = cell.y - 1;
            neigh_cell->x++;

            if (neigh_cell->x > cell.x + 1) {
                neigh_cell->x = cell.x - 1;
                return false;
            }
        }
    }

    return true;
}

/* 
 * cell_periodic_bound_cond: Check periodic boundary condition for cell. 
 *                           Cell must be inside simulation box, 
 *                           on the other hand, atom must be treated as it is
 *                           and must not be pulled back into the central
 *                           simation box.
 */
static void cell_periodic_bound_cond(cell_t cell_orig, cell_t *cell_new,
                                     double *atom_shift_x, 
                                     double *atom_shift_y, 
                                     double *atom_shift_z)
{
    /*
    cell_new->x = (cell_orig.x + cell_Lx) % cell_Lx;
    cell_new->y = (cell_orig.y + cell_Ly) % cell_Ly;
    cell_new->z = (cell_orig.z + cell_Lz) % cell_Lz;
    */

    /* Another way: */
    if (cell_orig.x < 0) {
        cell_new->x = cell_orig.x + cell_Lx;
        *atom_shift_x = -Lx;
    } else if (cell_orig.x >= cell_Lx) {
        cell_new->x = cell_orig.x - cell_Lx;
        *atom_shift_x = Lx;
    } else {
        cell_new->x = cell_orig.x;
        *atom_shift_x = 0;
    }

    if (cell_orig.y < 0) {
        cell_new->y = cell_orig.y + cell_Ly;
        *atom_shift_y = -Ly;
    } else if (cell_orig.y >= cell_Ly) {
        cell_new->y = cell_orig.y - cell_Ly;
        *atom_shift_y = Ly;
    } else {
        cell_new->y = cell_orig.y;
        *atom_shift_y = 0;
    }

    if (cell_orig.z < 0) {
        cell_new->z = cell_orig.z + cell_Lz;
        *atom_shift_y = -Lz;
    } else if (cell_orig.z >= cell_Lz) {
        cell_new->z = cell_orig.z - cell_Lz;
        *atom_shift_y = Lz;
    } else {
        cell_new->z = cell_orig.z;
        *atom_shift_z = 0;
    }
}

/* init_atom_in_cell: Init atom before scan in the cell. */
static bool init_atom_in_cell(cell_t cell, int *atom_i)
{
    if (head[cell.scal] != LIST_END) {
        *atom_i = head[cell.scal];
        return true;
    }
    else {
        *atom_i = -1;
        return false;
    }
}

/* scan_atom_in_cell: Scan in the cell: get list of atoms in the cell. */
static bool scan_atom_in_cell(cell_t cell, int *atom_i)
{
    if (list[*atom_i] != LIST_END) {
        *atom_i = list[*atom_i];
        return true;
    } else {
        *atom_i = -1;
        return false;
    }
}
