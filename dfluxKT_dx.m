function [ fx ] = dfluxKT_dx( ~, ~, ~ )
global n
    fx(:,1,1) = ones(n,1);
    fx(:,1,2) = zeros(n,1);
    fx(:,1,3) = zeros(n,1);
    
    fx(:,2,1) = zeros(n,1);
    fx(:,2,2) = ones(n,1);
    fx(:,2,3) = zeros(n,1);
    
    fx(:,3,1) = zeros(n,1);
    fx(:,3,2) = zeros(n,1);
    fx(:,3,3) = ones(n,1);
end