/*
 * tersoff.c: Tersoff potential.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "md.h"
#include "util.h"
#include "tersoff1.h"

/* 
 * Constants for atoms interactions (f_repulsive and f_attractive). 
 */

/* Parameters for silicon: */
#define A           3264.7      /* eV */
#define B           95.373      /* eV */
#define lambda1     3.2394      /* A^-1 */
#define lambda2     1.3258      /* A^-1 */
#define lambda3     lambda2
#define alpha       0.0
#define beta        0.33675
#define c           4.8381
#define d           2.0417
#define h           0.0
#define n           22.956
#define R           3.0         /* A */
#define D           0.2         /* A */
/* #define DOUBLE_CMP_EPS DBL_EPSILON */ /* Value from float.h. */
#define DOUBLE_CMP_EPS 1e-6 /* Value from float.h. */

static inline double t1_b(atom_t *ai, atom_t *aj, double rij);
static inline double t1_zeta(atom_t *ai, atom_t *aj, double rij);
static inline double t1_g(atom_t *ai, atom_t *aj, atom_t *ak);
static inline double t1_cos_theta(atom_t *ai, atom_t *aj, atom_t *ak);
static inline double t1_a(atom_t *ai, atom_t *aj, double rij);
static inline double t1_eta(atom_t *ai, atom_t *aj, double rij);
static inline double t1_f_cutoff(double rij);
static inline double t1_f_repulsive(double rij);
static inline double t1_f_attractive(double rij);
static inline double t1_distance(atom_t *ai, atom_t *aj);
static inline int    t1_iszero(double x);

/* tersoff1_forces: Computes forces. */
void tersoff1_forces()
{
    int i, j;
    double rij;           /* Distance between atom i and atom j */
    double Uij;           /* Bond energy between atom i and atom j */
    double U = 0;         /* Total energy. */
    double f_cutoff_val;
    
    for (i = 0; i < natoms; i++) {
        atom[i].ax = 0.0;
        atom[i].ay = 0.0;
        atom[i].az = 0.0;
    }

    for (i = 0; i < natoms - 1; i++) {
        for (j = i + 1; j < natoms; j++) {
            rij = t1_distance(&atom[i], &atom[j]);
            f_cutoff_val = t1_f_cutoff(rij);

            if (!t1_iszero(f_cutoff_val)) {
                /* If potential is enough large. */
                Uij = f_cutoff_val * (t1_a(&atom[i], &atom[j], rij) * 
                                      t1_f_repulsive(rij) + 
                                      t1_b(&atom[i], &atom[j], rij) * 
                                      t1_f_attractive(rij));
                U += Uij;
            }
        }
    }
    
    /*
    Andersen:
    Clean f[], a[]
    
    for i = 1 to n do
        epot(i) = 0
        
        for j = ii = 1 to neib(i) do
            if rij >= S then continue -- R_cutoff
            f_a =
            df_a
            f_r = 
            df_r =  
            fc() =
            dfc() = 
            
            // three-body term - b_ij
            for k = iii = 1 to neib(i)
                if k == j then cont.
                
                cosijk = 
                gz = 
                fc = 
                dfc = 
                dtau_j() = 
                dtau_k() = 
                dcsi_j() = 
                dcsi_k() = 
            end for    

            b_ij = ...
            a_ij = ...            
            
            epot(i) += ...
            
            f(i) +=  
            f(j) -=                
            
            a(i) += 
            a(j) -= 
            
            // three-body forces
            for k = iii = 1 to neib(i)
                f(i)
                f(j)
                f(k)
                a(i)
                a(j)
                a(k) 
            end for
        end for
    end for        
    
    */
    
}

/* */

/* t1_b: */
static inline double t1_b(atom_t *ai, atom_t *aj, double rij)
{
    return pow(1 + pow(beta, n) * 
           pow(t1_zeta(ai, aj, rij), n), -1/(2 * n));
}

/* t1_zeta: */
static inline double t1_zeta(atom_t *ai, atom_t *aj, double rij)
{
    double rik, f_cutoff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        atom_t *ak = &atom[k];

        if ((ak != ai) && (ak != aj)) {
            rik = t1_distance(ai, &atom[k]);
            f_cutoff_val = t1_f_cutoff(rik);

            if (f_cutoff_val > DOUBLE_CMP_EPS) {
                sum += f_cutoff_val * t1_g(ai, aj, &atom[k]) * 
                       exp(pow(lambda3, 3) * pow(rij - rik, 3));
            }
        }
    }
    return sum;
}

/* t1_g: */
double t1_g(atom_t *ai, atom_t *aj, atom_t *ak)
{
    return 1 + pow(c, 2) / pow(d, 2) -
           pow(c, 2) / (pow(d, 2) + (h - t1_cos_theta(ai, aj, ak)));
}

/* t1_cos_theta: */
inline double t1_cos_theta(atom_t *ai, atom_t *aj, atom_t *ak)
{
    return
    ((aj->x - ai->x) * (ak->x - ai->x) + 
     (aj->y - ai->y) * (ak->y - ai->y) +
     (aj->z - ai->z) * (ak->z - ai->z)) / 
    (sqrt(pow(aj->x - ai->x, 2) + pow(aj->y - ai->y, 2) + 
          pow(aj->z - ai->z, 2)) *
     sqrt(pow(ak->x - ai->x, 2) + pow(ak->y - ai->y, 2) + 
          pow(ak->z - ai->z, 2)));
}

/* t1_a: */
double t1_a(atom_t *ai, atom_t *aj, double rij)
{
    return pow(1 + pow(alpha, n) * 
                   pow(t1_eta(ai, aj, rij), n), -1 / (2 * n));
}

/* t1_eta: */
double t1_eta(atom_t *ai, atom_t *aj, double rij)
{
    double rik, f_cutoff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        atom_t *ak = &atom[k];

        if ((ak != ai) && (ak != aj)) {
            rik = t1_distance(ai, ak);
            f_cutoff_val = t1_f_cutoff(rik);

            if (f_cutoff_val > DOUBLE_CMP_EPS) {
                sum += f_cutoff_val * exp(pow(lambda3, 3) * 
                                          pow(rij - rik, 3));
            }
        }
    }
    return sum;
}

/* t1_f_cutoff: Smooth cutoff function. */
inline double t1_f_cutoff(double rij)
{
    if (rij < R - D) {
        return 1;
    } else if ((rij > R - D) && (rij < R + D)) {
        return 1 / 2 - 1 / 2 * sin(M_PI / 2 * (rij - R) / D);
    } else { /* if (rij > R + D) */
        return 0;
    }
}

/* t1_f_repulsive: Repulsive pair potential. */
inline double t1_f_repulsive(double rij)
{
    return A * exp(-lambda1 * rij);
}

/* t1_f_attractive: Attractive pair potential. */
inline double t1_f_attractive(double rij)
{
    return -B * exp(-lambda2 * rij);
}

/* t1_distance: Compute distance between atom i and atom j. */
static inline double t1_distance(atom_t *ai, atom_t *aj)
{
    return sqrt(pow(ai->x - aj->x, 2) + 
                pow(ai->y - aj->y, 2) + 
                pow(ai->z - aj->z, 2));
}

/* t1_iszero: Test floating point x is not a zero. */
static inline int t1_iszero(double x)
{
    return x > DOUBLE_CMP_EPS;
}
