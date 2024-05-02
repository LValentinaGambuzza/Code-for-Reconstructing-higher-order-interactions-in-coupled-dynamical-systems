%%% MATLAB code of the paper
%%% Reconstructing higher-order interactions in coupled dynamical systems
%%% By Federico Malizia, Alessandra Corso,Lucia Valentina Gambuzza
%%% Giovanni Russo, Vito Latora, Mattia Frasca
%%% Nature Communications -  NCOMMS-23-20846A

%%% This file finds a solution for the problem of reconstructing the
%%% interactions in the Lotka-Volterra model under the assumptions that
%%% the derivatives of the state variables are known
%%% It generates panel c of Figure 1

clear all

%parameters of the units
R =[0.6099;    0.6177;    0.8594;    0.8055;    0.5767;    0.1829; 0.2399];
K =[88.7647;    3.8387;   49.5002;   17.6248;   97.8894;   71.5568;   50.5467];
N=7;
A=[0 0 0 0 0 0 0;
    0 0 0 0 0.8 0 0;
    0 0.1 0 0 0 0 0;
    -0.4 0 0 0 0 -0.3 0;
    0 0 0 0.7 0 0 0;
    0 0 0 0 0 0 -0.7;
    0 -0.8 0 0.2 0 0 0];
A=A.*R./K;
g=digraph(A);
TriangleList=[2 3 7 0.0062;
    4 1 6 0.0016];

% initial conditions for the oscillators
x0=[30 45 32 50 55 30 40]';

% Maximum integration time
tmax = 20;

%initialization of the errors
err0 = [];
err1 = [];
err2 = [];
err4 = [];

%range of M/H
imin = 0.2;
imax = 4;
istep=0.05;

%number of unknowns for each node
H=(N-1)+(N-1)*(N-2)/2;

for iM = imin:istep:imax
    disp(iM)
    M = H*iM;
    dt = tmax/M;
    T = 0:dt:tmax;

    % Accurate integration of the equation
    options = odeset('abstol',1e-14,'reltol',1e-14);
    [T,X]=ode45(@(t,x) LVeqComplete(t,x,A,R,K,TriangleList),T,x0,options);

    nt = length(T);

    I1 = 2:nt-1; % Internal steps
    I2 = 3:nt-2; % Internal internal points

    Y1 = (X(I1+1,:)-X(I1,:))/dt;
    Y2 = (X(I1+1,:)-X(I1-1,:))/(2*dt);
    Y4 = (4*(X(I2+1,:)-X(I2-1,:))/(2*dt)-(X(I2+2,:)-X(I2-2,:))/(4*dt))/3;
    F = zeros(nt-2,N);

    for n=0:nt-1
        pop = X(n+1,:);
        F(n+1,:) = LVeqComplete(n*dt,pop',A,R,K,TriangleList);
    end
    Ftrue=F;
    
    A0=zeros(H,N); %reconstructed matrix
    AA=zeros(H,N); %"true" matrix
    TriangleList0=[2 3 7 0; 4 1 6 0];
    for n=0:nt-1
        pop = X(n+1,:);
        F(n+1,:) = LVeqComplete(n*dt,pop',zeros(N,N),R,K,TriangleList0);
    end

    % Zero order method
    Phi = zeros(nt-2,H);
    for i=1:N
        %pairwise terms j=1:N-1
         Phi(:,1:N-1) = X(I1,i).*X(I1,[1:i-1 i+1:N]);
         AA(1:N-1,i)=A(i,[1:i-1 i+1:N]);
        %higher-order terms j=N:H
        vtemp=zeros(nt-2,(N-1)*(N-2)/2);
        itemp=1;
        for ii1=[1:i-1 i+1:N]
            for jj1=[1:i-1 i+1:N]
                if jj1>ii1,
                    vtemp(:,itemp)=X(I1,i).*X(I1,ii1).*X(I1,jj1);
                    x1=find(TriangleList(:,1)==i);
                    x2=find(TriangleList(:,2)==ii1);
                    x3=find(TriangleList(:,3)==jj1);
                    if ((~isempty(x1))&&(~isempty(x2))&&(~isempty(x3)))
                        if ((x1==x2)&&(x2==x3)),
                            b=TriangleList(x1,4);
                            AA(N-1+itemp,i)=b;
                        end
                    end
                    itemp=itemp+1;
                end
            end
        end
        Phi(:,N:N-1+(N-1)*(N-2)/2)=vtemp;
        %solution
        Yi = Ftrue(I1,i)-F(I1,i);
        Zi=lsqminnorm(Phi,Yi);
        A0(:,i) = Zi;
    end
    err0 = [err0,norm(AA-A0,"fro")];
    
end

irange = imin:istep:imax;
normA = norm(AA,"fro");

figure,semilogy(irange,err0/normA);
xlabel('M/H')
ylabel('E')

figure,plot(irange,err0/normA);
xlabel('M/H')
ylabel('E')