    function yt = fct_pde3(t,y)
    global n alpha beta gamma v z0 zL D2 D1 z
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

     %% limiteur de pentes
    %Xz=KT_centered_limiter_fz(2,n,z,t,[C T],@fluxKT,@dfluxKT_dx);
    %Xz=KT_centered_limiter_fz_order1(2,n,z,t,[C T],@fluxKT,@dfluxKT_dx);

    %Cz = Xz(1:n)';
    %Tz = Xz(n+1:2*n)';
    
    %fz =KT_centered_limiter_fz(3,n,z,t,[u1, u2, u3],@fluxKT,@dfluxKT_dx);
    %fz =KT_centered_limiter_fz_order1(3,n,z,t,[u1, u2, u3],@fluxKT,@dfluxKT_dx);
    %u1z = fz(1:n)';
    %u2z = fz(n+1:2*n)';
    %u3z = fz(2*n+1:3*n)';

    %u1z=vanleer_slope_limiter_fz(z,n,t,u1,'flux','dflux_dx');
    %u2z=vanleer_slope_limiter_fz(z,n,t,u2,'flux','dflux_dx');
    %u3z=vanleer_slope_limiter_fz(z,n,t,u3,'flux','dflux_dx');

    %u1z=mc_slope_limiter_fz(z,n,t,u1,'flux','dflux_dx');
    %u2z=mc_slope_limiter_fz(z,n,t,u2,'flux','dflux_dx');
    %u3z=mc_slope_limiter_fz(z,n,t,u3,'flux','dflux_dx');
    
    %u1z=minmod_slope_limiter_fz(z,n,t,u1,'flux','dflux_dx');
    %u2z=minmod_slope_limiter_fz(z,n,t,u2,'flux','dflux_dx');
    %u3z=minmod_slope_limiter_fz(z,n,t,u3,'flux','dflux_dx');

    %u1z=smart_slope_limiter_fz(z,n,t,u1,'flux','dflux_dx');
    %u2z=smart_slope_limiter_fz(z,n,t,u2,'flux','dflux_dx');
    %u3z=smart_slope_limiter_fz(z,n,t,u3,'flux','dflux_dx');

    u1z=superbee_slope_limiter_fz(z,n,t,u1,'flux','dflux_dx');
    u2z=superbee_slope_limiter_fz(z,n,t,u2,'flux','dflux_dx');
    u3z=superbee_slope_limiter_fz(z,n,t,u3,'flux','dflux_dx');
   
%
%   Ecriture des equations
%

    u1t = -v*u1z+(1+alpha)*(1-u1).*u2.^2+beta*(1-u1).*u3.^2;
    u2t= -v*u2z+(1-alpha)*(1-u1).*u2.^2-gamma*u2;
    u3t= -v*u3z+beta*(1-u1).*u3.^2+2*alpha*(1-u1).*u2.^2-(gamma/beta)*u3;
    
    yt = [u1t; u2t; u3t];
    end
