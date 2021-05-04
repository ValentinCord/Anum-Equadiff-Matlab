    close all
    clear variables global
    
    global n alpha beta gamma v d1 d2 d3 z0 zL D1 z dz
    
% constantes 
    alpha = 0.065;
    beta = 2;
    gamma = 0.025;
    v = 1;
    d1 = 0.0001;
    d2 = 0.0002;
    d3 = 0.0002;
    
%
%   grille spatiale
%
    z0 = 0;
    zL = 1;
    n = 501;
    dz = (zL - z0)/(n-1);
    z = (z0:dz:zL)';
    
% conditions initiales
    ui1 = ones(n,1);
    ui2 = zeros(n,1);
    ui3 = zeros(n,1);
    
    y = [ui1; ui2; ui3];
    
%
%   matrices de diff�rentiation
%
   %D1 = two_point_upwind_D1(z,1);           %bien
   %D1 = three_point_centered_D1(z);          %Non
   %D1 = three_point_upwind_D1(z,1);          %moyen 
   %D1 = four_point_biased_upwind_D1(z,1);      %très bon 
   %D1 = four_point_upwind_D1(z,1);          %Non
   D1 = five_point_centered_D1(z);           
   %D1 = five_point_biased_upwind_D1(z,1);      $ pas mal
   
   %D2 n'est plus nécéssaire

%
%   instants de visualisation
%
    tI = 0;
    tF = 2;
    dt = 0.1;
    time = (tI:dt:tF);
    nt = length(time);
    
%   Option 
    
fx(:,1,1) = ones(n,1);
    fx(:,1,2) = zeros(n,1);
    fx(:,2,1) = zeros(n,1);
    fx(:,2,2) = ones(n,1);
    
    s2 = zeros(n);
    s2(1:length(s2) +1:numel(s2))=1;

    reltol = 1e-6;
    abstol = 1e-6;
    
    % sparse reduit la mémoire 
    % spones spones remplace les valleurs non nul par 1 

    jpattern(1:n,1:n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),2) + diag(ones(1,n-2),-2);
    jpattern(1:n,n+1:2*n) = eye(n);
    jpattern(1:n,2*n+1:3*n) = eye(n);
    
    jpattern(n+1:2*n,1:n) = eye(n);
    jpattern(n+1:2*n,n+1:2*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),2) + diag(ones(1,n-2),-2);
    jpattern(1:n,2*n+1:3*n) = eye(n);
    
    jpattern(2*n+1:3*n,1:n) = eye(n);
    jpattern(2*n+1:3*n,n+1:2*n) = eye(n);
    jpattern(2*n+1:3*n,2*n+1:3*n) = eye(n) + diag(ones(1,n-1),1) + diag(ones(1,n-1),-1) + diag(ones(1,n-2),2) + diag(ones(1,n-2),-2);
    
    %jpattern = sparse(spones(jpattern));
    options = odeset('RelTol', reltol, 'AbsTol', abstol, 'jpattern', jpattern);
    %options = odeset('RelTol', reltol, 'AbsTol', abstol);
%
%   int�gration temporelle
%
    tic
    %[timeout,yout]= ode45(@fct_pde,time,y, options);
    %[timeout,yout]= ode23(@fct_pde,time,y);
    [timeout,yout]= ode15s(@fct_pde3,time,y, options);
    %[timeout,yout]= ode113(@fct_pde,time,y, options);
    %[timeout,yout]= ode23s(@fct_pde,time,y, options);
    %[timeout,yout]= ode23t(@fct_pde,time,y, options);
    %[timeout,yout]= ode23tb(@fct_pde,time,y, options);
    
    toc
 %affichage
    u1 = yout(:, 1:n);
    u2 = yout(:,n+1:2*n);
    u3 = yout(:,2*n+1:3*n);
    
    for i =1:length(timeout)
      u1(i,1) = 0;
      u2(i,1) = 1;
      u3(i,1) = 0;
      u1(i,n) = u1(i,n-1);
      u2(i,n) = u2(i,n-1);
      u3(i,n) = u3(i,n-1);
    end
      
    titre = ['NbrPoint= ', num2str(n), ' Temps=  ', num2str(toc), 'sec'];
    subplot (3,1,1);
    plot(z,u1);
    xlabel ('z');
    ylabel('u1(z,t)');
    
    title(titre);
    subplot (3,1,2);
    plot(z,u2);
    xlabel ('z');
    ylabel('u2(z,t)');
    
    subplot (3,1,3);
    plot(z,u3);
    xlabel ('z');
    ylabel('u3(z,t)');
    
