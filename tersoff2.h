/*
 * tersoff2.h: Tersoff 2 potential.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#ifndef TERSOFF2_H
#define TERSOFF2_H

void tersoff2_energy();

inline double t2_f_cutoff(int i, int j, double rij);
inline double t2_f_repulsive(int i, int j, double rij);
inline double t2_f_attractive(int i, int j, double rij);

inline double t2_b(int i, int j);
inline double t2_zeta(int i, int j);
inline double t2_g(int i, int j, int k);
inline double t2_cos_theta(int i, int j, int k);

inline double t2_distance(int atom_i, int atom_j);
inline double t2_distance_noPBC(int atom_i, int atom_j);
inline int    iszero(double x);

#endif /* TERSOFF2_H */
