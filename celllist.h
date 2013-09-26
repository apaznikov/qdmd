/*
 * cellist.h: Linked cell list implementation
 * 
 * Copyright (C) 2013 Alexey Paznikov
 */

#ifndef CELLLIST_H
#define CELLLIST_H

typedef struct cell_s {
    int x;
    int y;
    int z;
    int scal;
} cell_t;

void linked_cell_init();
void linked_cell_finalize();

void linked_cell_update();
bool cell_update_test();

void init_cell(cell_t *cell);
bool scan_cells(cell_t *cell);
void init_neigh_cell(cell_t cell,  cell_t *neigh_cell);
bool scan_neigh_cells(cell_t cell, cell_t *neigh_cell);

bool init_atom_in_cell(cell_t cell, int *atom_i);
bool scan_atom_in_cell(cell_t cell, int *atom_i);

void cell_vec_to_scal(cell_t *cell);
void cell_periodic_bound_cond(cell_t cell_orig, cell_t *cell_new);

#endif /* end of include guard: CELLLIST_H */
