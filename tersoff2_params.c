/*
 * tersoff2_parameters.c: Parameters for tersoff 2 potential.
 * 
 * Copyright (C) 2013 Alexey Paznikov
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>

#include "md.h"
#include "tersoff2_params.h"
#include "util.h"

typedef struct param_single_s {
    double beta;
    double n;
    double h;
    double m;

    double mass;
} param_single_t;

typedef struct param_pair_s {
    double A;
    double B;
    double lambda;
    double mu;
    double R;
    double S;
    double chi;
    double epsi;
} param_pair_t;

typedef struct param_3bound_s {
    double sigma;
    double c;
    double d;
    double omega;
    double delta;
} param_3bound_t;

enum {
    NSPECIES = 4
};

param_single_t carbon_single_param, 
               silicon_single_param,
               germanium_single_param,
               hydrogen_single_param;

param_pair_t   carbon_pair_param, 
               silicon_pair_param,
               germanium_pair_param,
               hydrogen_pair_param;

param_3bound_t carbon_3bound_param,
               silicon_3bound_param,
               germanium_3bound_param,
               hydrogen_3bound_param;

param_single_t param_all_single[NSPECIES];
param_pair_t   param_all_pair[NSPECIES];
param_pair_t   **param_pair_precomp;
param_3bound_t param_all_3bound[NSPECIES];

#define chi_Si_C  0.97760e0
#define chi_Si_Ge 1.00061e0
#define chi_Si_H  1.0485e0

static void init_single_param();
static void init_pair_precomp_param();

static double lambda_ij_precomp(int spec_i, int spec_j);
static double lambda_precomp(int spec_i);
static double mu_ij_precomp(int spec_i, int spec_j);
static double mu_precomp(int spec_i);
static double Aij_precomp(int spec_i, int spec_j);
static double A_precomp(int spec_i);
static double Bij_precomp(int spec_i, int spec_j);
static double B_precomp(int spec_i);
static double Rij_precomp(int spec_i, int spec_j);
static double R_precomp(int spec_i);
static double Sij_precomp(int spec_i, int spec_j);
static double S_precomp(int spec_i);
static double chi_ij_precomp(int spec_i, int spec_j);

/* tersoff2_init_param: Init parameters for atoms. */
void tersoff2_init_param()
{
    init_single_param();
    init_pair_precomp_param();
}

/* init_single_param: Init single-element parameters. */
static void init_single_param()
{
    /* Single */

    carbon_single_param.beta     = 1.5724e-7;
    carbon_single_param.n        = 7.2751e-1;
    carbon_single_param.h        = -5.7058e-1;
    carbon_single_param.m        = 1.0e0;
    carbon_single_param.mass     = 12.011e0;

    silicon_single_param.beta    = 1.1000e-6;
    silicon_single_param.n       = 7.8734e-1;
    silicon_single_param.h       = -5.9825e-1;
    silicon_single_param.m       = 1.0e0;
    silicon_single_param.mass    = 28.0855e0;

    germanium_single_param.beta  = 9.0166e-7;
    germanium_single_param.n     = 7.5627e-1;
    germanium_single_param.h     = -4.3884e-1;
    germanium_single_param.m     = 1.0e0;
    germanium_single_param.mass  = 72.61e0;

    hydrogen_single_param.beta   = 1.0e0;
    hydrogen_single_param.n      = 1.e0;
    hydrogen_single_param.h      = -1.0e0;
    hydrogen_single_param.m      = 1.6094e0;
    hydrogen_single_param.mass   = 1.0079e0;

    /* Pair */

    carbon_pair_param.A          = 1.3936e3;
    carbon_pair_param.B          = 3.4670e2;
    carbon_pair_param.lambda     = 3.4879e0;
    carbon_pair_param.mu         = 2.2119e0;
    carbon_pair_param.R          = 1.8e0;
    carbon_pair_param.S          = 2.1e0;
    carbon_pair_param.chi        = 1.0e0;  
    carbon_pair_param.epsi       = 1.0e0;

    silicon_pair_param.A         = 1.8308e3;
    silicon_pair_param.B         = 4.7118e2;
    silicon_pair_param.lambda    = 2.4799e0;
    silicon_pair_param.mu        = 1.7322e0;
    silicon_pair_param.R         = 2.7e0;
    silicon_pair_param.S         = 3.0e0; 
    silicon_pair_param.chi       = 1.0e0;    
    silicon_pair_param.epsi      = 1.0e0;

    germanium_pair_param.A       = 1.7690e3;
    germanium_pair_param.B       = 4.1923e2;
    germanium_pair_param.lambda  = 2.4451e0;
    germanium_pair_param.mu      = 1.7047e0;
    germanium_pair_param.R       = 2.8e0;
    germanium_pair_param.S       = 3.1e0;  
    germanium_pair_param.chi     = 1.0e0;      
    germanium_pair_param.epsi    = 1.0e0;

    hydrogen_pair_param.A        = 8.007e1;
    hydrogen_pair_param.B        = 3.138e1;
    hydrogen_pair_param.lambda   = 4.2075e0;
    hydrogen_pair_param.mu       = 1.7956e0;
    hydrogen_pair_param.R        = 1.1e0;
    hydrogen_pair_param.S        = 1.7e0;    
    hydrogen_pair_param.chi      = 1.0e0;        
    hydrogen_pair_param.epsi     = 1.0e0;

    /* 3-bound */

    carbon_3bound_param.sigma    = 0.0e0;
    carbon_3bound_param.c        = 3.8049e4;
    carbon_3bound_param.d        = 4.3840e0;
    carbon_3bound_param.omega    = 1.e0;
    carbon_3bound_param.delta    = 0.e0;     

    silicon_3bound_param.sigma   = 0.0e0;
    silicon_3bound_param.c       = 1.0039e5;
    silicon_3bound_param.d       = 1.6217e1;
    silicon_3bound_param.omega   = 1.e0;
    silicon_3bound_param.delta   = 0.e0;

    germanium_3bound_param.sigma = 0.0e0;
    germanium_3bound_param.c     = 1.0643e5;
    germanium_3bound_param.d     = 1.5652e1;
    germanium_3bound_param.omega = 1.e0;
    germanium_3bound_param.delta = 0.e0;

    hydrogen_3bound_param.sigma  = 3.0e0;
    hydrogen_3bound_param.c      = 0.e0;
    hydrogen_3bound_param.d      = 1.e0;
    hydrogen_3bound_param.omega  = 4.e0;
    hydrogen_3bound_param.delta  = 0.e0;

    /* Set correspondence between species and species numbers. */

    param_all_single[0] = carbon_single_param;
    param_all_single[1] = silicon_single_param;
    param_all_single[2] = germanium_single_param;
    param_all_single[3] = hydrogen_single_param;

    param_all_pair[0]   = carbon_pair_param;
    param_all_pair[1]   = silicon_pair_param;
    param_all_pair[2]   = germanium_pair_param;
    param_all_pair[3]   = hydrogen_pair_param;

    param_all_3bound[0] = carbon_3bound_param;
    param_all_3bound[1] = silicon_3bound_param;
    param_all_3bound[2] = germanium_3bound_param;
    param_all_3bound[3] = hydrogen_3bound_param;
}

/* init_pair_precomp_param: Precompute pair parameters. */
static void init_pair_precomp_param()
{
    int spec_i, spec_j;

    param_pair_precomp = xmalloc(sizeof(param_pair_t*) * natomspecies);
    for (spec_i = 0; spec_i < natomspecies; spec_i++) {
        param_pair_precomp[spec_i] = xmalloc(sizeof(param_pair_t)*natomspecies);
    }

    for (spec_i = 0; spec_i < natomspecies; spec_i++) {
        for (spec_j = 0; spec_j < natomspecies; spec_j++) {
            param_pair_precomp[spec_i][spec_j].lambda = 
                lambda_ij_precomp(spec_i, spec_j);
            param_pair_precomp[spec_i][spec_j].mu = 
                mu_ij_precomp(spec_i, spec_j);
            param_pair_precomp[spec_i][spec_j].A = 
                Aij_precomp(spec_i, spec_j);
            param_pair_precomp[spec_i][spec_j].B = 
                Bij_precomp(spec_i, spec_j);
            param_pair_precomp[spec_i][spec_j].R = 
                Rij_precomp(spec_i, spec_j);
            param_pair_precomp[spec_i][spec_j].S = 
                Sij_precomp(spec_i, spec_j);
            param_pair_precomp[spec_i][spec_j].chi = 
                chi_ij_precomp(spec_i, spec_j);
        }
    }
}

/*
 * Pair parameters.
 */

/* lambda_ij: */
inline double lambda_ij(int i, int j)
{
    return param_pair_precomp[atomspecies[i]][atomspecies[j]].lambda;
}

/* lambda_ij_precomp: */
static double lambda_ij_precomp(int spec_i, int spec_j)
{
    return (lambda_precomp(spec_i) + lambda_precomp(spec_j)) / 2;
}

/* lambda_precomp: */
static double lambda_precomp(int spec_i)
{
    return param_all_pair[spec_i].lambda;
}

/* mu_ij: */
inline double mu_ij(int i, int j)
{
    return param_pair_precomp[atomspecies[i]][atomspecies[j]].mu;
}

/* mu_ij_precomp */
static double mu_ij_precomp(int spec_i, int spec_j)
{
    return (mu_precomp(spec_i) + mu_precomp(spec_j)) / 2;
}

/* mu_precomp */
static double mu_precomp(int spec_i)
{
    return param_all_pair[spec_i].mu;
}

/* Aij: */
inline double Aij(int i, int j)
{
    return param_pair_precomp[atomspecies[i]][atomspecies[j]].A;
}

/* Aij_precomp: */
static double Aij_precomp(int spec_i, int spec_j)
{
    return sqrt(A_precomp(spec_i) * A_precomp(spec_j));
}

/* A_precomp: */
static double A_precomp(int spec_i)
{
    return param_all_pair[spec_i].A;
}

/* Bij: */
inline double Bij(int i, int j)
{
    return param_pair_precomp[atomspecies[i]][atomspecies[j]].B;
}

/* Bij_precomp: */
static double Bij_precomp(int spec_i, int spec_j)
{
    return sqrt(B_precomp(spec_i) * B_precomp(spec_j));
}

/* B_precomp: */
static double B_precomp(int spec_i)
{
    return param_all_pair[spec_i].B;
}

/* Rij: */
inline double Rij(int i, int j)
{
    return param_pair_precomp[atomspecies[i]][atomspecies[j]].R;
}

/* Rij_precomp: */
static double Rij_precomp(int spec_i, int spec_j)
{
    return sqrt(R_precomp(spec_i) * R_precomp(spec_j));
}

/* R_precomp: */
static double R_precomp(int spec_i)
{
    return param_all_pair[spec_i].R;
}

/* Sij: */
inline double Sij(int i, int j)
{
    return param_pair_precomp[atomspecies[i]][atomspecies[j]].S;
}

/* Sij_precomp: */
static double Sij_precomp(int spec_i, int spec_j)
{
    return sqrt(S_precomp(spec_i) * S_precomp(spec_j));
}

/* S_precomp: */
static double S_precomp(int spec_i)
{
    return param_all_pair[spec_i].S;
}

/* chi: */
inline double chi(int i, int j)
{
    if (i == j) {
        return 1.0;
    } else {
        return param_pair_precomp[atomspecies[i]][atomspecies[j]].chi;
    }
}

/* chi_ij_precomp: line 122 in tersoff.f90 */
static double chi_ij_precomp(int spec_i, int spec_j)
{
    if (strcmp(atomspecies_names[spec_i], "H") && 
        strcmp(atomspecies_names[spec_j], "H")) {

        if (!strcmp(atomspecies_names[spec_i], "C") || 
            !strcmp(atomspecies_names[spec_j], "C")) {
            return chi_Si_C;

        } else if (!strcmp(atomspecies_names[spec_i], "GE") || 
                   !strcmp(atomspecies_names[spec_j], "GE")) {
            return chi_Si_Ge;
        }
    }

    /* NOTE strange logic in tersoff.f90 */
    return chi_Si_H;
}

/*
 * Single parameters
 */

/* omega: Parameter omega is available in the future. Here it equals 1. */
inline double omega_ij(int i, int j)
{
    return 1.0;
}

/* c: */
inline double c(int i)
{
    return param_all_3bound[atomspecies[i]].c;
}

/* d: */
inline double d(int i)
{
    return param_all_3bound[atomspecies[i]].d;
}

/* h: */
inline double h(int i)
{
    return param_all_single[atomspecies[i]].h;
}

/* n: */
inline double n(int i)
{
    return param_all_single[atomspecies[i]].n;
}

/* beta: */
inline double beta(int i)
{
    return param_all_single[atomspecies[i]].beta;
}

/* mass: */
inline double mass(int i)
{
    return param_all_single[atomspecies[i]].mass;
}

