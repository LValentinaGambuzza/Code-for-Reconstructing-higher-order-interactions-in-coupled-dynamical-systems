%%% MATLAB code of the paper
%%% Reconstructing higher-order interactions in coupled dynamical systems
%%% By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
%%% Giovanni Russo, Vito Latora, Mattia Frasca
%%% Nature Communications -  NCOMMS-23-20846A

%%% This file finds a solution for the problem of reconstructing the
%%% interactions in the Roessler model under the assumptions that
%%% the derivatives of the state variables are known
%%% It generates the results of panel d of Figure 2

clear all

%remove the following three lines, if not using parallel computing 
if isempty(gcp('nocreate'))
     parpool('local',16);
end

load('ZackaryNet.mat');

for i45=1:length(EdgeList)
    A(EdgeList(i45,1),EdgeList(i45,2))=1;
    A(EdgeList(i45,2),EdgeList(i45,1))=1;
end

TriangleList=closedtriangles;

%parameters of the coupling
k=1e-4;
kD=1e-5;

% initial conditions for the oscillators
xoold =(30*rand(N,1)-15)/5;
yoold =(30*rand(N,1)-15)/5;
zoold =(40*rand(N,1)-5)/5;
x0=[xoold; yoold; zoold];

% parameters of integration
dt=0.01;
tmax=100;
T=0:dt:tmax;
options=odeset('abstol',1e-12,'reltol',1e-12);

%range of M/H
imin = 100;
imax = 1000;
istep=5;
irange=round(logspace(1,3,16));

iMvalori=irange;

%initialization of errors
nvaloriM=length(irange)
err0 = zeros(nvaloriM,1);
err00 = zeros(nvaloriM,1);
err2 = zeros(nvaloriM,1);
err3 = zeros(nvaloriM,1);


%number of unknowns for each node
H=(N-1)+(N-1)*(N-2)/2;

%change parfor into for, if not using parallel computing 
parfor iMindice = 1:length(iMvalori)
    iM=iMvalori(iMindice);
    disp(iM)
    M = H*iM;
    dt = tmax/M;
    T = 0:dt:tmax;

    % Accurate integration of the equation
    [T,X]=ode45(@(t,x) roessler_hoi(t,x,EdgeList,TriangleList),T,x0,options);
    
    nt = length(T);

    I1 = 2:nt-1; % Internal steps
    I2 = 3:nt-2; % Internal internal points

    F = zeros(nt-2,3*N);

    for n=0:nt-1
        pop = X(n+1,:);
        F(n+1,:)=roessler_hoi(n*dt,pop',EdgeList,TriangleList);
    end
    Ftrue=F;
    
    Y4 = (4*(X(I2+1,:)-X(I2-1,:))/(2*dt)-(X(I2+2,:)-X(I2-2,:))/(4*dt))/3;

    
    A4=zeros(H,N); %reconstructed matrix with OLS, fourth order approx
    A5=zeros(H,N); %reconstructed matrix with NNLS, fourth order approx
    AA=zeros(H,N); %vera
    EdgeList0=[1 1; 2 2];
    TriangleList0=[1 1 1; 2 2 2];
    for n=0:nt-1
        pop = X(n+1,:);
        F(n+1,:) = roessler_hoi(n*dt,pop',EdgeList0,TriangleList0);
    end

   
    % Fourth order method
    Phi = zeros(nt-4,H);
    for i=1:N
        %pairwise terms j=1:N-1
        Phi(:,1:N-1) = k*(X(I2,[1:i-1 i+1:N])-X(I2,i));
        AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
        %h.o.i. terms j=N:H
        vtemp=zeros(nt-4,(N-1)*(N-2)/2);
        itemp=1;
        for ii1=[1:i-1 i+1:N]
            for jj1=[1:i-1 i+1:N]
                if jj1>ii1,
                    vtemp(:,itemp)=kD*(X(I2,ii1).*X(I2,jj1).^2+X(I2,jj1).*X(I2,ii1).^2-2*X(I2,i).^3);
                    RowIdx = find(ismember(TriangleList, [i ii1 jj1],'rows'));
                    if (~isempty(RowIdx))
                        AA(N-1+itemp,i)=1;
                    end
                    RowIdx = find(ismember(TriangleList, [ii1 i jj1],'rows'));
                    if (~isempty(RowIdx))
                        AA(N-1+itemp,i)=1;
                    end
                    RowIdx = find(ismember(TriangleList, [ii1 jj1 i],'rows'));
                    if (~isempty(RowIdx))
                        AA(N-1+itemp,i)=1;
                    end
                    itemp=itemp+1;
                end
            end
        end
        Phi(:,N:N-1+(N-1)*(N-2)/2)=vtemp;
        
        
        %solution
        Yi = Y4(:,i)-F(I2,i);
        Zi4=lsqminnorm(Phi,Yi,1e-12);
        
        Zi5 = lsqnonneg(Phi,Yi);
        
        A4(:,i) = Zi4;
        A5(:,i) = Zi5;
        
    end
    
    
    normA = norm(AA,"fro");
    err4(iMindice,1)=norm(AA-A4,"fro")/normA;
    err5(iMindice,1)=norm(AA-A5,"fro")/normA;
    
    
end

save results.mat

figure,loglog(irange,err4,irange,err5);
xlabel('M/H')
ylabel('E')