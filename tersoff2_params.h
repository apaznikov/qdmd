/*
 * tersoff2_parameters.h: Parameters for tersoff 2 potential.
 * 
 * Copyright (C) 2013 Alexey Paznikov
 */

#ifndef TERSOFF2_PARAMETERS_H
#define TERSOFF2_PARAMETERS_H

void tersoff2_param_init();
void tersoff2_param_finalize();

inline double lambda_ij(atom_t *ai, atom_t *aj);
inline double mu_ij(atom_t *ai, atom_t *aj);
inline double Aij(atom_t *ai, atom_t *aj);
inline double Bij(atom_t *ai, atom_t *aj);
inline double Rij(atom_t *ai, atom_t *aj);
inline double Sij(atom_t *ai, atom_t *aj);
inline double omega_ij(atom_t *ai, atom_t *aj);
inline double c(atom_t *ai);
inline double d(atom_t *ai);
inline double h(atom_t *ai);
inline double n(atom_t *ai);
inline double chi(atom_t *ai, atom_t *aj);
inline double beta(atom_t *ai);
inline double mass(atom_t *ai);

#endif /* POTENTIAL_TERSOFF2_H */
