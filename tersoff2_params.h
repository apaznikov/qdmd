/*
 * tersoff2_parameters.h: Parameters for tersoff 2 potential.
 * 
 * Copyright (C) 2013 Alexey Paznikov
 */

#ifndef TERSOFF2_PARAMETERS_H
#define TERSOFF2_PARAMETERS_H

void tersoff2_param_init();
void tersoff2_param_finalize();

inline double lambda_ij(int i, int j);
inline double mu_ij(int i, int j);
inline double Aij(int i, int j);
inline double Bij(int i, int j);
inline double Rij(int i, int j);
inline double Sij(int i, int j);
inline double omega_ij(int i, int j);
inline double c(int i);
inline double d(int i);
inline double h(int i);
inline double n(int i);
inline double chi(int i, int j);
inline double beta(int i);
inline double mass(int i);

#endif /* POTENTIAL_TERSOFF2_H */
