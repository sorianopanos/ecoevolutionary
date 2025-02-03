#!/bin/bash
#$ -t 1-1
#$ -tc 500

cc MCMF_evolution.c -o evolution -lm
mu="0.143"
interacciones="10"
nrealizations="200"
alpha0="0.00"
path="../Results/trajectories_stoch_crossimmunity/"
Dinfectivity="0.000"
flag_model="stochastic"
flag_execution="trajectories"
mutation_probability=".10"
cross_immunity="1"
R0="3"
Dimmunity="0.03"
humanostotal="10000"
Niter="1500"

./evolution ${mu} ${Dinfectivity} ${Dimmunity} ${interacciones} ${R0} ${nrealizations} ${alpha0} ${mutation_probability} ${cross_immunity} ${humanostotal} ${Niter} ${flag_model} ${flag_execution} ${path}