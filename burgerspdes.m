    function xt = burgerspdes(t,u)
    global mu
    global n z0 zL D2 D1
    t;
%
%   conditions aux limites
%
    u(1) = burgers_exact(z0,t);
    u(n) = burgers_exact(zL,t);
%
%   dérivées temporelles :
%
    xt = - u.*(D1*u) + mu*D2*u;
    end
