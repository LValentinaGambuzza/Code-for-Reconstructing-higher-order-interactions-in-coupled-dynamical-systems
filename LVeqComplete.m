function dNdt=LVeqComplete(t,Npop,A,R,K,TriangleList)
%%% MATLAB code of the paper
%%% Reconstructing higher-order interactions in coupled dynamical systems
%%% By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
%%% Giovanni Russo, Vito Latora, Mattia Frasca
%%% Nature Communications -  NCOMMS-23-20846A

%%% This file contains the equations of the Lotka-Volterra model
%%% including pairwise and higher-order interactions



%pairwise coupling
coup_rete=A*Npop;

%h.o.i. coupling
coup_simplicial=zeros(7,1);
[mtrianglelist,ntrianglelist]=size(TriangleList);
for ii=1:mtrianglelist
    i1=TriangleList(ii,1);
    i2=TriangleList(ii,2);
    i3=TriangleList(ii,3);
    b=TriangleList(ii,4);
    coup_simplicial(i1)=coup_simplicial(i1)+b*Npop(i2,1)*Npop(i3,1);
end

dNdt=Npop.*R.*(1-Npop./K)+Npop.*coup_rete+Npop.*coup_simplicial;        
