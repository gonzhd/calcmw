# mwcalc

Calculation of the Moleular Weight with corresponding uncertainty 
for a given Molecular Formula. 
The molecular weight is calculated by summing over the atomic weights of the 
constituents elements. 

    $MW=\sum(N_i*atw_i)$
    
The atomic weights, in turn, are calculated as the weighted average of the atomic
masses using the isotopic distribution as weights.

    $atw_i=\sum(w_j*atm_j)$
    
Uncertainties are calculated by uncertainty propagation of these two equations.