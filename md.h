/*
 * main.h: Main module.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#ifndef MAIN_H
#define MAIN_H

enum {
    ATOMTYPE_NAME_SIZE = 2
};

typedef struct {
    /* Atoms positions, velocities and accelerations */
    double x, y, z;
    double vx, vy, vz;
    double ax, ay, az;
} atom_t;

extern int natoms;
extern atom_t *atom;

extern int *atomspecies;
extern int natomspecies;
extern char **atomspecies_names;

/* Experimental box size. */
extern double Lx, Ly, Lz;

/* Timesteps */
int ntimesteps;
int timestep;
double timestep_dt;

#ifndef _STDBOOL_H
typedef enum 
{
    true  = 1,
    false = 0
} bool;
#endif /* _STDBOOL_H */

void md_initialize(const char *cls_file, const char *input_file);
void md_finalize();
void md_run();

#endif /* MAIN_H */
