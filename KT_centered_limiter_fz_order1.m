      function [fz] = KT_centered_limiter_fz_order1(ne,n,z,t,x,flux,dflux_dx)
%...
%...  The MatMol Group (2009)
%...
%...  function KT_centered_slope_limiter_order1 is a order-one version of
%...  KT_centered_limiter : see KT_centered_limiter for explanations
%...
%...  argument list
%...
%...     ne         number of scalar equations (input)
%...
%...     n          number of grid points in the z domain including the
%...                boundary points (input)
%...
%...     z          independent variable (input) : z(n) which dimensions
%...
%...                must be 1 x n
%...
%...     t          time (input)
%...
%...     x          dependent variable (input) : x(n,ne)
%...
%...     flux       matlab function which computes f(x) : f(n,ne)
%...                call: f = fun(t,x) where flux = 'fun'     (input)
%...
%...     dflux_dx   matlab function which computes the derivative of f with
%...                respect to x : dfdx(n,ne,ne)
%...                call: f_x = fun(t,x) where dflux_dx = 'fun'
%...
%...
%... addendum by Ph Saucez (June 2014) : a new variable called zloc is
%... introduced : zloc contains the elements of z and is allways a row
%... vector, whatever is z : a row or a column vector

    zloc(1:n) = z(1:n);

      for i=1:ne
          xp(1:n-1,i) = x(2:n,i);
          xm(1:n-1,i) = x(1:n-1,i);
      end
      
      dfludxp = feval(dflux_dx,ne,t,xp);
      dfludxm = feval(dflux_dx,ne,t,xm);

      
      for i=1:n-1
          ajac(1:ne,1:ne) = dfludxp(i,1:ne,1:ne);
          bajac           = balance(ajac);
          vp(1:ne)        = eig(bajac);
          ajac(1:ne,1:ne) = dfludxm(i,1:ne,1:ne);
          bajac           = balance(ajac);
          vp(ne+1:2*ne)   = eig(bajac);
          aspeed(i,1)     = max(abs(vp));
      end

      flzp = feval(flux,ne,t,xp);
      flzm = feval(flux,ne,t,xm);
      flz = feval(flux,ne,t,x);

      for j=1:ne
          fz(1,j)= (flz(2,j)-flz(1,j))/(zloc(2)-zloc(1));
          fz(2:n-1,j)=(flzp(2:n-1,j)-flzm(1:n-2,j)...
                       -aspeed(2:n-1,1).*(xp(2:n-1,j)-xm(2:n-1,j))...
                       +aspeed(1:n-2,1).*(xp(1:n-2,j)-xm(1:n-2,j)))./(zloc(3:n)-zloc(1:n-2))';
          fz(n,j)= (flz(n,j)-flz(n-1,j))/(zloc(n)-zloc(n-1));
      end
