    function f = fct(t,u1,u2,u3)
    %global mu
    global n z0 zL D2 D1 alpha beta gamma v d1 d2 d3
    t;
%
%   conditions aux limites
%
    u1(:,1)=0;
    u2(:,1)=1;
    u3(:,1)=0;
    %D1*u1(:,n)=0;
    %D1*u2(:,n)=0;
    %D1*u3(:,n)=0;
%
%   dérivées temporelles :
%
    
    f = [d1*D2*u1' - v*D1*u1' + (1+alpha)*(1-u1)*u2^2 + beta*(1-u1)*u3^2;
    d2*D2*u2' - v*D1*u2' +(1-alpha)*(1-u1)*u2^2-gamma*u2;
    d3*D2*u3'-v*D1*u3'+beta*(1-u1)*u3^2+2*alpha*(1-u1)*u2^2-gamma*u3/beta];
    end