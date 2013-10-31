/*
 * tersoff2.h: Tersoff 2 potential.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov, Alexey Paznikov
 */

#ifndef TERSOFF2_H
#define TERSOFF2_H

void tersoff2_energy();

inline double t2_f_cutoff(atom_t *ai, atom_t *aj, double rij);
inline double t2_f_repulsive(atom_t *ai, atom_t *aj, double rij);
inline double t2_f_attractive(atom_t *ai, atom_t *aj, double rij);

inline double t2_b(atom_t *ai, atom_t *aj);
inline double t2_zeta(atom_t *ai, atom_t *aj);
inline double t2_g(atom_t *ai, atom_t *aj, atom_t *ak);
inline double t2_cos_theta(atom_t *ai, atom_t *aj, atom_t *ak);

inline double t2_distance(atom_t *ai, atom_t *aj);
inline double t2_distance_noPBC(atom_t *ai, atom_t *aj);
inline int    iszero(double x);

#endif /* TERSOFF2_H */
