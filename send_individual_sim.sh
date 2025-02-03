#!/bin/bash
#$ -t 1-1
#$ -tc 500

cc MCMF_evolution.c -o evolution -lm
mu="0.143"  ### set the recovery rate
k="10" ## set number of contacts k
nrealizations="200" ### number of realizations to average
alpha0="0.00" ### put a different value from 0 to activate the recovery-transmissibility tradeoff
path="../Results/trajectories_stoch_crossimmunity/" ### set the path to store the results
Dinfectivity="0.000" ### changes in transmissibility
flag_model="stochastic" ### stochastic (normal model), deterministic (model with determinsitic changes in evolution and transmissibility), tradeoff (model incorporating the transmissibility-recovery tradeoff)
flag_execution="trajectories" ## trajectories for the time evolution of epidemic curves, endemic for the fraction of simulations surviving after Niter time steps 
mutation_probability="1" #### probability that the variant of an infected individual mutates at one time step
cross_immunity="1" #### Always set to 1
R0="3" ## basic reproduction number of the wild-type variant
Dimmunity="0.03" ## changes in antigenic position 
population="10000" ### size of the population 
Niter="1500" ## number of days simulated in each epidemic trajectory

./evolution ${mu} ${Dinfectivity} ${Dimmunity} ${k} ${R0} ${nrealizations} ${alpha0} ${mutation_probability} ${cross_immunity} ${population} ${Niter} ${flag_model} ${flag_execution} ${path}
