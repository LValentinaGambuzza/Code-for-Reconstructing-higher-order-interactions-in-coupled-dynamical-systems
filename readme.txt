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
including pairwise and higher-order interactions.
It generates panel b of Figure 1

LVexact.m
This file finds a solution for the problem of reconstructing the
interactions in the Lotka-Volterra model under the assumptions that
the derivatives of the state variables are known.
It generates panel c of Figure 1

LVapproximations.m
This file finds a solution for the problem of reconstructing the
interactions in the Lotka-Volterra model using different approximations
for the derivatives
It generates panel d of Figure 1

roessler_hoi.m
This file contains the equations of the Roessler model
including pairwise and higher-order interactions

ZackaryNet.mat
MATLAB file containing the data to simulate the Zackary club higher-order structure

signal_lasso.m
This file is the MATLAB implementation of the signal lasso algorithm
based on the paper L. Shi, C. Shen, L. Jin, Q. Shi, Z. Wang, 
and S. Boccaletti, Physical Review Research 3, 043210 (2021)

RosslerExact.m
This file finds a solution for the problem of reconstructing the
interactions in the Roessler model under the assumptions that
the derivatives of the state variables are known.
It generates panel b of Figure 2

RosslerApprox.m
This file finds a solution for the problem of reconstructing the
interactions in the Roessler model using different approximations
for the derivatives.
It generates the results of panel d of Figure 2