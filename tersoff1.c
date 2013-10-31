/*
 * tersoff.c: Tersoff potential.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov
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

static inline double t1_b(int i, int j, double rij);
static inline double t1_zeta(int i, int j, double rij);
static inline double t1_g(int i, int j, int k);
static inline double t1_cos_theta(int i, int j, int k);
static inline double t1_a(int i, int j, double rij);
static inline double t1_eta(int i, int j, double rij);
static inline double t1_f_cutoff(double rij);
static inline double t1_f_repulsive(double rij);
static inline double t1_f_attractive(double rij);
static inline double t1_distance(int i, int j);
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
        ax[i] = 0.0;
        ay[i] = 0.0;
        az[i] = 0.0;
    }

    for (i = 0; i < natoms - 1; i++) {
        for (j = i + 1; j < natoms; j++) {
            rij = t1_distance(i, j);
            f_cutoff_val = t1_f_cutoff(rij);

            if (!t1_iszero(f_cutoff_val)) {
                /* If potential is enough large. */
                Uij = f_cutoff_val * (t1_a(i, j, rij) * t1_f_repulsive(rij) + 
                                      t1_b(i, j, rij) * t1_f_attractive(rij));
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
inline double t1_b(int i, int j, double rij)
{
    return pow(1 + pow(beta, n) * 
           pow(t1_zeta(i, j, rij), n), -1/(2 * n));
}

/* t1_zeta: */
inline double t1_zeta(int i, int j, double rij)
{
    double rik, f_cutoff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        if ((k != i) && (k != j)) {
            rik = t1_distance(i, k);
            f_cutoff_val = t1_f_cutoff(rik);

            if (f_cutoff_val > DOUBLE_CMP_EPS) {
                sum += f_cutoff_val * t1_g(i, j, k) * exp(pow(lambda3, 3) * 
                                                      pow(rij - rik, 3));
            }
        }
    }
    return sum;
}

/* t1_g: */
inline double t1_g(int i, int j, int k)
{
    return 1 + pow(c, 2) / pow(d, 2) -
           pow(c, 2) / (pow(d, 2) + (h - t1_cos_theta(i, j, k)));
}

/* t1_cos_theta: */
inline double t1_cos_theta(int i, int j, int k)
{
    return
    ((x[j] - x[i]) * (x[k] - x[i]) + 
     (y[j] - y[i]) * (y[k] - y[i]) +
     (z[j] - z[i]) * (z[k] - z[i])) / 
    (sqrt(pow(x[j] - x[i], 2) + pow(y[j] - y[i], 2) + pow(z[j] - z[i], 2)) *
     sqrt(pow(x[k] - x[i], 2) + pow(y[k] - y[i], 2) + pow(z[k] - z[i], 2)));
}

/* t1_a: */
inline double t1_a(int i, int j, double rij)
{
    return pow(1 + pow(alpha, n) * 
                   pow(t1_eta(i, j, rij), n), -1 / (2 * n));
}

/* t1_eta: */
inline double t1_eta(int i, int j, double rij)
{
    double rik, f_cutoff_val, sum = 0;
    int k;

    for (k = 0; k < natoms; k++) {
        if ((k != i) && (k != j)) {
            rik = t1_distance(i, k);
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
static inline double t1_distance(int i, int j)
{
    return sqrt(pow(x[i] - x[j], 2) + 
                pow(y[i] - y[j], 2) + 
                pow(z[i] - z[j], 2));
}

/* t1_iszero: Test floating point x is not a zero. */
static inline int t1_iszero(double x)
{
    return x > DOUBLE_CMP_EPS;
}
