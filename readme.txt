MATLAB code of the paper
Reconstructing higher-order interactions in coupled dynamical systems
By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
Giovanni Russo, Vito Latora, Mattia Frasca
Nature Communications -  NCOMMS-23-20846A

LVeqComplete.m
This file contains the equations of the Lotka-Volterra model
including pairwise and higher-order interactions

LVmodel.m
This file integrates the equations of the Lotka-Volterra model
including pairwise and higher-order interactions
It generates panel b of Figure 1

LVexact.m
This file finds a solution for the problem of reconstructing the
interactions in the Lotka-Volterra model under the assumptions that
the derivatives of the state variables are known
It generates panel c of Figure 1


LVapproximations.m
This file finds a solution for the problem of reconstructing the
interactions in the Lotka-Volterra model using different approximations
for the derivatives
It generates panel d of Figure 1