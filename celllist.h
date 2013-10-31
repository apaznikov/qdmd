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

void celllist_scan_atoms(void (*func) (int atom_i, int atom_j));

void linked_cell_init();
void linked_cell_finalize();

void linked_cell_update();
bool cell_update_test();


#endif /* end of include guard: CELLLIST_H */
