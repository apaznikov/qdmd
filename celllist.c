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

/* linked_cell_init: Init head and all constants. */
void linked_cell_init()
{
    const double cellsize = r_cutoff + delta_r;

    int cell;

    /* Number of cells by coordinates */
    cell_Lx = floor(Lx / cellsize);
    cell_Ly = floor(Ly / cellsize);
    cell_Lz = floor(Lz / cellsize);

    ncells = cell_Lx * cell_Ly * cell_Lz;

    /* temp */
    cell_LyLz = cell_Ly * cell_Lz;   

    /* Size of a one cell. */
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
void cell_vec_to_scal(cell_t *cell)
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
void init_cell(cell_t *cell)
{
    cell->x = 0;
    cell->y = 0;
    cell->z = -1; /* -1 for first increment in scan_cells */
}

/* scan_cells: Scan over cells. Goes to the next cell. */
bool scan_cells(cell_t *cell)
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
void init_neigh_cell(cell_t cell, cell_t *neigh_cell)
{
    neigh_cell->x = cell.x - 1;
    neigh_cell->y = cell.y - 1;
    neigh_cell->z = cell.z - 1 - 1; /* -1 for first increment */
}

/* scan_neigh_cells: */
bool scan_neigh_cells(cell_t cell, cell_t *neigh_cell)
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

/* cell_periodic_bound_cond: Check periodic boundary condition for cell. */
void cell_periodic_bound_cond(cell_t cell_orig, cell_t *cell_new)
{
    if (cell_orig.x < 0) {
        cell_new->x = cell_orig.x + cell_Lx;
    } else if (cell_orig.x >= cell_Lx) {
        cell_new->x = cell_orig.x - cell_Lx;
    } else {
        cell_new->x = cell_orig.x;
    }

    if (cell_orig.y < 0) {
        cell_new->y = cell_orig.y + cell_Ly;
    } else if (cell_orig.y >= cell_Ly) {
        cell_new->y = cell_orig.y - cell_Ly;
    } else {
        cell_new->y = cell_orig.y;
    }

    if (cell_orig.z < 0) {
        cell_new->z = cell_orig.z + cell_Lz;
    } else if (cell_orig.z >= cell_Lz) {
        cell_new->z = cell_orig.z - cell_Lz;
    } else {
        cell_new->z = cell_orig.z;
    }
}

/* init_atom_in_cell: Init atom before scan in the cell. */
bool init_atom_in_cell(cell_t cell, int *atom_i)
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
bool scan_atom_in_cell(cell_t cell, int *atom_i)
{
    if (list[*atom_i] != LIST_END) {
        *atom_i = list[*atom_i];
        return true;
    } else {
        *atom_i = -1;
        return false;
    }
}
