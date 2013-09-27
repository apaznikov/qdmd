/*
 * md.c: Molecular dynamics main module.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "md.h"
#include "tersoff1.h"
#include "tersoff2.h"
#include "tersoff2_forces.h"
#include "tersoff2_params.h"
#include "celllist.h"
#include "util.h"

enum {
    INTPUTFILE_LINESIZE_MAX  = 1024,
    INTPUTFILE_PARAMSIZE_MAX = 1024
};

int natoms;
atom_t *atom;

int natomspecies;
int *atomspecies;
char **atomspecies_names;

double Lx, Ly, Lz;

/* Timesteps */
int ntimesteps = 500;
int timestep;
double timestep_dt = 5E-15;

static void allocate_memory();
static void free_memory();

static void read_atoms(const char *filename);
static void read_inputfile(const char *input_file);
static void read_inputfile_species(FILE *fin, const char *line);
static void read_inputfile_thermobox(FILE *fin);

static void integrate();
static void compute_pos();
static void compute_velocities();

/* 
 * md_run: Runs molecular dynamic simulation. 
 */
void md_run()
{
    double t;

    printf("Run modelling\n");

    tersoff2_param_init();

    linked_cell_init();
    linked_cell_update();
    
    /* Initialize simulation parameters */
    /* cutoff ... */
    
    /* Evaluate properties: E, preasure, ... */
    
    /* Compute forces for r0, v0, a0 */

    /* FIXME: it's for debug! uncomment! */
    tersoff2_energy();
    tersoff2_forces();
    
    for (timestep = 1; timestep <= ntimesteps; timestep++) {
        t = timestep * timestep_dt;

        printf("step %6d timestep: %E sec.\n", timestep, t);
        integrate();

        if (cell_update_test()) {
            linked_cell_update();
        }
    }

    tersoff2_param_finalize();
    linked_cell_finalize();
}

/* 
 * md_initialize: Initializes MD module: loads initial atoms positions, 
 *                velocities, setups parameters of potentials. 
 * Parameters:
 *   cls_file - file in CLSMAN format (initial positions, 
 *              velocities, atom species)
 *   input_file - file with modelling parameters
 */
void md_initialize(const char *cls_file, const char *input_file)
{
    read_atoms(cls_file);
    read_inputfile(input_file);
}

/* md_finalize: Finalizes MD module. */
void md_finalize()
{
    free_memory();
}

/* allocate_memory: Allocates memory for global arrays. */
static void allocate_memory()
{
    atom = xmalloc(sizeof(atom_t) * natoms);
    atomspecies = xmalloc(sizeof(*atomspecies) * natoms);
}

/* free_memory: */
static void free_memory()
{
    int species;

    free(atom);
    free(atomspecies);

    for (species = 0; species < natomspecies; species++) {
        free(atomspecies_names[species]);
    }
    free(atomspecies_names);
}

/* 
 * read_atoms: Reads atoms from file (CLSMAN format)
 *
 * Format:
 *   natoms NUM? NUM?
 *   line?
 *   periodic cell
 *   x y z atom_species vx vy vz
 *   ...
 * 
 */
static void read_atoms(const char *filename)
{
    enum {
        LINESIZE_MAX = 1024
    };
    char line[LINESIZE_MAX];
    int lineno, i;
    FILE *fin;
    
    printf("Read atoms from file '%s'\n", filename);

    if ( (fin = fopen(filename, "r")) == NULL)
        exit_error("Can't open file '%s'", filename);
    
    /* Read natoms */
    if (fgets(line, LINESIZE_MAX, fin) == NULL)
        exit_error("Incorrect format of '%s'", filename);
    line[strlen(line) - 1] = '\0';
  
    if (sscanf(line, "%d", &natoms) < 1)
        exit_error("Can't read number of atoms in '%s'", filename);

    /* Read title */
    if (fgets(line, LINESIZE_MAX, fin) == NULL)
        exit_error("Incorrect format of '%s'", filename);

    /* Read periodic info */
    if (fgets(line, LINESIZE_MAX, fin) == NULL)
        exit_error("Incorrect format of '%s'", filename);
    
    /* Read atoms */
    allocate_memory(natoms);
    lineno = 3;
    i = 0;
    while (fgets(line, LINESIZE_MAX, fin) != NULL) {
        lineno++;
        /* Read: x  y  z  atom_species  vx  vy  vz */
        if (sscanf(line, "%lf %lf %lf %d %lf %lf %lf", 
                   &atom[i].x, &atom[i].y, &atom[i].z, &atomspecies[i], 
                   &atom[i].vx, &atom[i].vy, &atom[i].vz) < 7)
        {                                                 
            exit_error("Can't read atom %d info (line %d)", i, lineno);
        }
        i++;
    }
    fclose(fin);
    
    printf("Loaded %d atoms\n", i);    
    if (i != natoms)
        exit_error("Number of loaded atoms is not equal to natoms (line 1)");
}

/* read_inputfile: Read input file. */
static void read_inputfile(const char *input_file)
{
    FILE *fin;
    char line[INTPUTFILE_LINESIZE_MAX];
    char param[INTPUTFILE_PARAMSIZE_MAX];
    
    printf("Read input from file '%s'\n", input_file);

    if ( (fin = fopen(input_file, "r")) == NULL)
        exit_error("Can't open file '%s'", input_file);

    while (fgets(line, INTPUTFILE_LINESIZE_MAX, fin) != NULL) {

        sscanf(line, "%s", param);

        if (!strcmp(param, "SPECIES")) {
            read_inputfile_species(fin, line);
        } else if (!strcmp(param, "THERMO")) {
            sscanf(line + strlen("THERMO"), "%s", param);

            if (!strcmp(param, "BOX")) {
                read_inputfile_thermobox(fin);
            }
        }
    }

    fclose(fin);
}

/* read_inputfile_species: Read atom species from input file. */
static void read_inputfile_species(FILE *fin, const char *line)
{
    int species;
    char ch;

    fseek(fin, -strlen(line) + strlen("SPECIES"), SEEK_CUR);
    fscanf(fin, "%d", &natomspecies);
    atomspecies_names = xmalloc(sizeof(char *) * natomspecies);

    for (species = 0; species < natomspecies; species++) {
        atomspecies_names[species] = xmalloc(sizeof(char) * 
                                     ATOMTYPE_NAME_SIZE);
    }

    for (species = 0; species < natomspecies; species++) {
        fscanf(fin, "%s", atomspecies_names[species]);
    }

    do {
        ch = fgetc(fin);
    } while (ch != '\n');
}

/* read_inputfile_thermobox: Read sizes of experimental box. */
static void read_inputfile_thermobox(FILE *fin)
{
    double x1, y1, z1, x2, y2, z2;
    char line[INTPUTFILE_LINESIZE_MAX];

    if (fgets(line, INTPUTFILE_LINESIZE_MAX, fin) != NULL) {
        sscanf(line, "%lf%lf%lf%lf%lf%lf", &x1, &y1, &z1, &x2, &y2, &z2);
    }
    
    Lx = fabs(x1 - x2);
    Ly = fabs(y1 - y2);
    Lz = fabs(z1 - z2);
}

/* 
 * integrate: Implements Velocity-Verlet time integration algorithms. 
 */
static void integrate()
{
    compute_pos();
    
    tersoff2_forces();
    tersoff2_energy();
    
    compute_velocities();
}

/*
 * compute_pos: Current time is t and we known previous state: 
 *              r(t - dt), v(t - dt), a(t - dt)
 *              Compute current position based on Velocity Verlet algorithm
 *              and velocity v(t + 1/2 * dt)
 */
static void compute_pos()
{
    int ai;

    for (ai = 0; ai < natoms; ai++) {
        /* Compute position. */
        atom[ai].x = atom[ai].x + atom[ai].vx * timestep_dt + 
                     0.5 * atom[ai].ax * timestep_dt * timestep_dt; 
        atom[ai].y = atom[ai].y + atom[ai].vy * timestep_dt + 
                     0.5 * atom[ai].ay * timestep_dt * timestep_dt; 
        atom[ai].z = atom[ai].z + atom[ai].vz * timestep_dt + 
                     0.5 * atom[ai].az * timestep_dt * timestep_dt; 

        /* Periodic boundary conditions. */
        /* TODO: make function for PBC */
        if (atom[ai].x >= Lx / 2) {
            atom[ai].x = atom[ai].x - Lx;
        } else if (atom[ai].x < -Lx / 2) {
            atom[ai].x = atom[ai].x + Lx;
        }

        if (atom[ai].y >= Ly / 2) {
            atom[ai].y = atom[ai].y - Ly;
        } else if (atom[ai].y < -Ly / 2) {
            atom[ai].y = atom[ai].y + Ly;
        }

        if (atom[ai].z >= Lz / 2) {
            atom[ai].z = atom[ai].z - Lz;
        } else if (atom[ai].z < -Lz / 2) {
            atom[ai].z = atom[ai].z + Lz;
        }

        /* Update velocities based on a(t - dt) */
        atom[ai].vx = atom[ai].vx + 0.5 * atom[ai].ax * timestep_dt;
        atom[ai].vy = atom[ai].vy + 0.5 * atom[ai].ay * timestep_dt;
        atom[ai].vz = atom[ai].vz + 0.5 * atom[ai].az * timestep_dt;
    }
}

/* compute_velocities: Compute velocities v(t) based on a(t) */
static void compute_velocities()
{
    int ai;

    for (ai = 0; ai < natoms; ai++) {
        atom[ai].vx = atom[ai].vx + 0.5 * atom[ai].ax * timestep_dt;
        atom[ai].vy = atom[ai].vy + 0.5 * atom[ai].ay * timestep_dt;
        atom[ai].vz = atom[ai].vz + 0.5 * atom[ai].az * timestep_dt;
    }
}
