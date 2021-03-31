      function [fz] = KT_centered_limiter_fz(ne,n,z,t,x,flux,dflux_dx)
%...
%...  The MatMol Group (2009)
%...
%...  function KT_centered_limiter returns the first derivative, fz,
%...  of a flux-function vector f(x) over the spatial domain z0 < z < zL.
%...
%...  This function implements a centered scheme for convection-diffusion PDEs
%...  proposed by Kurganov and Tadmor (2000), which can be used as a general
%...  finite-volume method, independently of the eigenstructure of the problem.
%...  The only required information is an estimation of the velocity of the wave.
%...
%...  argument list
%...
%...     ne         number of equations (input)
%...
%...     n          number of grid points in the z domain including the
%...                boundary points (input)
%...
%...     z          independent variable (input) : z(n)
%...
%...     t          time (input)
%...
%...     x          dependent variable (input) : x(n,ne) = [x1 x2 ... xne]
%...
%...     flux       matlab function which computes f(x)= [f1(x1,...,xne) f2(x1,...,xne) ... fne(x1,...,xne)]
%...                The syntax is :
%...
%...                        function out = flux(ne,t,x)
%...                        n = length(x)
%...                        out(1:n,1) = f1(x1,...,xne);
%...                        ...
%...                        out(1:n,ne) = fne(x1,...,xne);
%...
%...                the variable t in the calling list makes possible to
%...                implement the case where the flux function f depends
%...                explicitly of t : f(t,x) in the place of f(x)
%...
%...     dflux_dx   matlab function which computes the derivative of f with
%...                respect to x : the output is a three-dimensionnal
%...                matrix out(1:n,1:ne,1:ne). The syntax is:
%...
%...                        function out = dflux_dx(ne,t,x)
%...                        n = length(x)
%...                        out(1:n,1,1)   = df1/dx1
%...                        ...
%...                        out(1:n,1,ne)  = df1/dxne
%...                        out(1:n,2,1)   = df2/dx1
%...                        ...
%...                        out(1:n,2,ne)  = df2/dxne
%...                        ...
%...                        ...
%...                        out(1:n,ne,1)  = dfne/dx1
%...                        ...
%...                        out(1:n,ne,ne) = dfne/dxne
%...                        
%...                the variable t in the calling list makes possible to
%...                implement the case where the dflux_dx function depends
%...                explicitly of t.
%...

%...  In compact form, for ne = 1, and on a uniform dz spaced gird,
%...  the resulting central scheme of Kurganov and Tadmor is
%...
%...               +           -             +           -
%...  df       [f(x     ) + f(x     )] - [f(x     ) + f(x     )] 
%...    j          j+1/2       j+1/2         j-1/2       j-1/2
%...  --- =    -------------------------------------------------
%...  dz                        2*dz
%...
%...             1            +       -                 +       -
%...         - ---- [a     *(x     - x     ) - a     *(x     - x     )]
%...           2*dz   j+1/2   j+1/2   j+1/2     j-1/2   j-1/2   j-1/2
%...
%... The following code generalizes that scheme to a system of ne pdes on a
%... non uniform grid.
%...
%...
%... The first loop constructs for each of the ne dependent variables
%...
%...     +               dz 
%...    x      = x    -  -- x        in xp 
%...     j+1/2    j+1     2  (j+1)z       
%...
%... and
%...
%...     -             dz 
%...    x      = x  +  -- x    in xm 
%...     j+1/2    j     2  jz       
%...
%... where
%...                  x - x      x   - x        x   - x
%...                   j   j-1    j+1   j-1      j+1   j
%...    x  = minmod[2*-------- , ---------- , 2*--------]
%...     jz           z - z      z   - z        z   - z
%...                   j   j-1    j+1   j-1      j+1   j
%...
%...                    
%... and minmod(a1,a2,a3,...) = min(aj) if all aj > 0
%...                             j   
%...                          = max(aj) if all aj < 0
%...                             j
%...                          = 0 otherwise
%...
%...
%... addendum by Ph Saucez (June 2014) : a new variable called zloc is
%... introduced : zloc contains the elements of z and is allways a row
%... vector, whatever is z : a row or a column vector

    zloc(1:n) = z(1:n);

    for i=1:ne
          xtmp(1:n) = x(1:n,i);
          
          xtmpz(1) = (xtmp(2)-xtmp(1))/(zloc(2)-zloc(1));

          xtmpz1(1:n-2) = 2*(xtmp(3:n)-xtmp(2:n-1))./(zloc(3:n)-zloc(2:n-1));
          xtmpz2(1:n-2) = (xtmp(3:n)-xtmp(1:n-2))./(zloc(3:n)-zloc(1:n-2));          
          xtmpz3(1:n-2) = 2*(xtmp(2:n-1)-xtmp(1:n-2))./(zloc(2:n-1)-zloc(1:n-2));
          testsgn       = sign(xtmpz1)+sign(xtmpz2)+sign(xtmpz3);
          for j=1:n-2
              if (testsgn(j) == 3) || (testsgn(j) == -3)
                  xtmpz(j+1)=sign(testsgn(j))*min([abs(xtmpz1(j)) abs(xtmpz2(j)) abs(xtmpz3(j))]);
              else
                  xtmpz(j+1)=0;
              end
          end

          xtmpz(n)    = (xtmp(n)-xtmp(n-1))/(zloc(n)-zloc(n-1));
          
          xtp(1:n-1)  = xtmp(2:n)-xtmpz(2:n).*(zloc(2:n)-zloc(1:n-1))/2;
          xtm(1:n-1)  = xtmp(1:n-1)+xtmpz(1:n-1).*(zloc(2:n)-zloc(1:n-1))/2;
          
          xp(1:n-1,i) = xtp(1:n-1);
          xm(1:n-1,i) = xtm(1:n-1);
      end
%...
%... The second loop computes the a      parameters, called local speed of propagation
%...                               j+1/2
%...
%...                                +                     -
%... with a     = max[rho(dflux_dx(x     )),rho(dflux_dx(x     ))]
%...       j+1/2                    j+1/2                 j+1/2 
%...
%... and where rho(A) is the spectral radius of A
%...
%...
      dfludxp = feval(dflux_dx,ne,t,xp);
      dfludxm = feval(dflux_dx,ne,t,xm);
      
      for i=1:n-1
          ajac(1:ne,1:ne) = dfludxp(i,1:ne,1:ne);
%...
%... balance is a diagonal scaling to improve eigenvalues accuracy
%...
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
%...      
%...                         df
%...                           j
%... The third loop computes --- : for j = 1 and j = n, simplified formulas are  
%...                         dz
%... required because central schemes are not practicables.
%...
      for j=1:ne
          fz(1,j)= (flz(2,j)-flz(1,j))/(zloc(2)-zloc(1));
          fz(2:n-1,j)=(flzp(2:n-1,j)-flzp(1:n-2,j)+flzm(2:n-1,j)-flzm(1:n-2,j)...
                       -aspeed(2:n-1,1).*(xp(2:n-1,j)-xm(2:n-1,j))...
                       +aspeed(1:n-2,1).*(xp(1:n-2,j)-xm(1:n-2,j)))./(zloc(3:n)-zloc(1:n-2))';
          fz(n,j)= (flz(n,j)-flz(n-1,j))/(zloc(n)-zloc(n-1));
      end
          