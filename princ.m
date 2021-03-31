    close all
    clear variables global
%    global mu
    global n z0 zL D2 D1 alpha beta gamma v d1 d2 d3
    
%
%   grille spatiale
%
    z0 = 0;
    zL = 1;
    n = 500;
    dz = (zL - z0)/(n-1);
    z = (z0:dz:zL)';
    v = 1;
    alpha = 0.065;
    beta = 2;
    d1 = 0.0001;
    d2 = 0.0002;
    d3 = 0.0002;
%
%   matrices de différentiation
%
   D1 = two_point_upwind_D1(z,v);
   D2 = three_point_centered_D2(z);
%
%   constantes du problème
%
%    mu = 0.005;
%
%   conditions initiales
%
    u1 = zeros(1,n);
    u2 = zeros(1,n);
    u3 = zeros(1,n);
%    u1(1,:)=1;
%    u2(1,:)=0;
%    u3(1,:)=0;
%
%   instants de visualisation
%
    dt = 0.1;
    time = (0:dt:1);
    nt = length(time);
%
%   intégration temporelle
%
    [timeout,yout]= ode45(@fct,time,[u1 u2 u3]);
    
    figure
    hold
    plot(z,yout,'.-k')
    
