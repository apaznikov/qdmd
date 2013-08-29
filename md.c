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
#include "util.h"

enum {
    INTPUTFILE_LINESIZE_MAX  = 1024,
    INTPUTFILE_PARAMSIZE_MAX = 1024
};


/* Atoms positions, velocities and accelerations */
int natoms;
double *x, *y, *z;
double *vx, *vy, *vz;
double *ax, *ay, *az;
int *atomspecies;

int natomspecies;
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

    tersoff2_init_param();
    
    /* Initialize simulation parameters */
    /* cutoff ... */
    
    /* Evaluate properties: E, preasure, ... */
    
    /* Compute forces for r0, v0, a0 */

    /* FIXME: it's for debug! uncomment!
     * tersoff2_energy(); */
    tersoff2_forces();
    
    for (timestep = 1; timestep <= ntimesteps; timestep++) {
        t = timestep * timestep_dt;
        printf("step %6d timestep: %E sec.\n", timestep, t);
        integrate();
    }
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
    x = xmalloc(sizeof(*x) * natoms);
    y = xmalloc(sizeof(*y) * natoms);
    z = xmalloc(sizeof(*z) * natoms);
    vx = xmalloc(sizeof(*vx) * natoms);
    vy = xmalloc(sizeof(*vy) * natoms);
    vz = xmalloc(sizeof(*vz) * natoms);
    ax = xmalloc(sizeof(*ax) * natoms);
    ay = xmalloc(sizeof(*ay) * natoms);
    az = xmalloc(sizeof(*az) * natoms);
    atomspecies = xmalloc(sizeof(*atomspecies) * natoms);
}

/* free_memory: */
static void free_memory()
{
    int species;

    free(x);
    free(y);
    free(z);
    free(vx);
    free(vy);
    free(vz);
    free(ax);
    free(ay);
    free(az);
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
                   &x[i], &y[i], &z[i], &atomspecies[i], 
                   &vx[i], &vy[i], &vz[i]) < 7)
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

/* read_inputfile_thermobox: */
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
    int i;

    for (i = 0; i < natoms; i++) {
        /* Compute position. */
        x[i] = x[i] + vx[i] * timestep_dt + 
               0.5 * ax[i] * timestep_dt * timestep_dt; 
        y[i] = y[i] + vy[i] * timestep_dt + 
               0.5 * ay[i] * timestep_dt * timestep_dt; 
        z[i] = z[i] + vz[i] * timestep_dt + 
               0.5 * az[i] * timestep_dt * timestep_dt; 

        /* Periodic boundary conditions. */
        if (x[i] >= Lx / 2) {
            x[i] = x[i] - Lx;
        } else if (x[i] < -Lx / 2) {
            x[i] = x[i] + Lx;
        }

        if (y[i] >= Ly / 2) {
            y[i] = y[i] - Ly;
        } else if (y[i] < -Ly / 2) {
            y[i] = y[i] + Ly;
        }

        if (z[i] >= Lz / 2) {
            z[i] = z[i] - Lz;
        } else if (z[i] < -Lz / 2) {
            z[i] = z[i] + Lz;
        }

        /* Update velocities based on a(t - dt) */
        vx[i] = vx[i] + 0.5 * ax[i] * timestep_dt;
        vy[i] = vy[i] + 0.5 * ay[i] * timestep_dt;
        vz[i] = vz[i] + 0.5 * az[i] * timestep_dt;
    }
}

/* compute_velocities: Compute velocities v(t) based on a(t) */
static void compute_velocities()
{
    int i;

    for (i = 0; i < natoms; i++) {
        vx[i] = vx[i] + 0.5 * ax[i] * timestep_dt;
        vy[i] = vy[i] + 0.5 * ay[i] * timestep_dt;
        vz[i] = vz[i] + 0.5 * az[i] * timestep_dt;
    }
}
