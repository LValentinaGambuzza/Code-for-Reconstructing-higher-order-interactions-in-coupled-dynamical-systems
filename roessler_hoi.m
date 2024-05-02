function dxdt=roessler_hoi(t,x,EdgeList,TriangleList)
%%% MATLAB code of the paper
%%% Reconstructing higher-order interactions in coupled dynamical systems
%%% By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
%%% Giovanni Russo, Vito Latora, Mattia Frasca
%%% Nature Communications -  NCOMMS-23-20846A

%%% This file contains the equations of the Roessler model
%%% including pairwise and higher-order interactions

[m1,n1]=size(x);

xold=x(1:m1/3,1);
yold=x((m1/3+1):m1/3*2,1);
zold=x((m1/3*2+1):end,1);

N=m1/3;

% parameters of the units
ar=0.2;
br=0.2;
cr=9;

%parameters of the coupling
k=1e-4;
kD=1e-5;

coup_rete=zeros(N,1);
coup_simplicial=zeros(N,1);

%pairwise terms
for ii=1:length(EdgeList)
    i1=EdgeList(ii,1);
    i2=EdgeList(ii,2);
    coup_rete(i1)=coup_rete(i1)+xold(i2)-xold(i1);
    coup_rete(i2)=coup_rete(i2)+xold(i1)-xold(i2);
end

%higher-order interactions
[mtrianglelist,ntrianglelist]=size(TriangleList);
for ii=1:mtrianglelist
    i1=TriangleList(ii,1);
    i2=TriangleList(ii,2);
    i3=TriangleList(ii,3);
    coup_simplicial(i1)=coup_simplicial(i1)+xold(i2)^2*xold(i3)-xold(i1)^3+xold(i2)*xold(i3)^2-xold(i1)^3;
    coup_simplicial(i2)=coup_simplicial(i2)+xold(i1)^2*xold(i3)-xold(i2)^3+xold(i1)*xold(i3)^2-xold(i2)^3;
    coup_simplicial(i3)=coup_simplicial(i3)+xold(i1)^2*xold(i2)-xold(i3)^3+xold(i1)*xold(i2)^2-xold(i3)^3;
end

dxdt1=-yold-zold+k*coup_rete+kD*coup_simplicial;
dydt1=xold+ar*yold;
dzdt1=br+zold.*(xold-cr);

dxdt=[dxdt1; dydt1; dzdt1];