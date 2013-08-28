/*
 * tersoff.c: Tersoff potential.
 * 
 * Copyright (C) 2013 Mikhail Kurnosov
 */

#include <stdio.h>
#include <stdlib.h>

#include "md.h"
#include "util.h"

/* forces_tersoff: Computes forces. */
void forces_tersoff()
{
    int i;
    
    for (i = 0; i < natoms; i++) {
        ax[i] = 0.0;
        ay[i] = 0.0;
        az[i] = 0.0;
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
