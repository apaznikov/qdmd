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
bool scan_cells(cell_t *cell);
bool scan_neigh_cells(cell_t cell, cell_t *neigh_cell);

#endif /* end of include guard: CELLLIST_H */
