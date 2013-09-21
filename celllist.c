/*
 * cellist.c: Linked cell list implementation
 * 
 * Copyright (C) 2013 Alexey Paznikov
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "md.h"
#include "util.h"

enum { END = -1 }; /* End of a list flag. */

static int *list, *head;

void linked_cell_init();
void linked_cell_finalize();

/* linked_cell_init: */
void linked_cell_init()
{
    /* TODO: how to choose cutoff? */
    /* Larger than this disance interactions can be neglected. */
    /* r_cutoff = 2.5: [Rapaport, 48], [Allen, Tildesley, 148]  */
    double r_cutoff = 2.5;        

    /* Number of cells by each coordinate. */
    double cell_Lx = floor(Lx / r_cutoff);
    double cell_Ly = floor(Ly / r_cutoff);
    double cell_Lz = floor(Lz / r_cutoff);

    double cell_LyLz = cell_Ly * cell_Lz;   /* temp */
    int ncells = cell_Lx * cell_Ly * cell_Lz;

    /* Size of a one cell. */
    double cell_rx = Lx / cell_Lx;
    double cell_ry = Ly / cell_Ly;
    double cell_rz = Lz / cell_Lz;

    int cell_coord_x, cell_coord_y, cell_coord_z, cell_i;
    int cell, atom;

    list = xmalloc(sizeof(*list) * natoms);
    head = xmalloc(sizeof(*head) * ncells);

    /* Reset headers to emtpy (end) */
    for (cell = 0; cell < ncells; cell++) {
        head[cell] = END;
    }

    /* Construct head and list */
    for (atom = 0; atom < natoms; atom++) {
        cell_coord_x = x[atom] / cell_rx;
        cell_coord_y = y[atom] / cell_ry;
        cell_coord_z = z[atom] / cell_rz;

        cell_i = cell_coord_x * cell_LyLz + 
                 cell_coord_y * cell_Ly + 
                 cell_coord_z;

        list[atom] = head[cell_i];
        head[cell_i] = atom;
    }
}

/* linked_cell_finalize(): */
void linked_cell_finalize()
{
    free(list);
    free(head);
}
