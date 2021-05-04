    function yt = fct_pde2(t,y)
    global n alpha beta gamma v D1
    t
    
% recuperation du vecteur y des équations u1, u2 et u3
    u1 = y(1:n);
    u2 = y(n+1:2*n);
    u3 = y(2*n+1:3*n);
    
%
%   conditions aux limites
%
% Dirichelet
    u1(1) = 0;
    u2(1) = 1;
    u3(1) = 0;
    
%Neumann
    u1(n) = u1(n-1);
    u2(n) = u2(n-1);
    u3(n) = u3(n-1);
%
%   Ecriture des equations
%   
    u1t = -v*D1*u1+(1+alpha)*(1-u1).*u2.^2+beta*(1-u1).*u3.^2;
    u2t= -v*D1*u2+(1-alpha)*(1-u1).*u2.^2-gamma*u2;
    u3t= -v*D1*u3+beta*(1-u1).*u3.^2+2*alpha*(1-u1).*u2.^2-(gamma/beta)*u3;
    
    yt = [u1t; u2t; u3t];
    end
