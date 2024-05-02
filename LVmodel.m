%%% MATLAB code of the paper
%%% Reconstructing higher-order interactions in coupled dynamical systems
%%% By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
%%% Giovanni Russo, Vito Latora, Mattia Frasca
%%% Nature Communications -  NCOMMS-23-20846A

%%% This file integrates the equations of the Lotka-Volterra model
%%% including pairwise and higher-order interactions
%%% This file generates panel b of Figure 1

clear all

%parameters of the units
 R =[0.6099;    0.6177;    0.8594;    0.8055;    0.5767;    0.1829; 0.2399];
 K =[88.7647;    3.8387;   49.5002;   17.6248;   97.8894;   71.5568;   50.5467];

N=7;
%pairwise terms
A=[0 0 0 0 0 0 0;
    0 0 0 0 0.8 0 0;
    0 0.1 0 0 0 0 0;
    -0.4 0 0 0 0 -0.3 0;
    0 0 0 0.7 0 0 0;
    0 0 0 0 0 0 -0.7;
    0 -0.8 0 0.2 0 0 0];
A=A.*R./K;

g=digraph(A);

%h.o.i. terms
TriangleList=[2 3 7 0.0062;
    4 1 6 0.0016];

% initial conditions for the oscillators
x0=[30 45 32 50 55 30 40]';

% Integration of the equations
tmax=40;
dt=0.01;
T=0:dt:tmax;
options=odeset('abstol',1e-12,'reltol',1e-12);
[T,X]=ode45(@(t,x) LVeqComplete(t,x,A,R,K,TriangleList),T,x0,options);
figure,plot(T,X),ylim([0 250])
xlabel('time')
ylabel('x_i(t)')
numberofdeadspecies=sum(abs(X(end,:))<0.01)