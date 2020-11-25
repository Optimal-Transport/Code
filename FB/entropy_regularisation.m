function [gam,obj,temp] = Entropy_Regularisation(m,n,c,eps)
%%
% Entropy_Regularisation(m,n,c,eps) : Executes the Entropy Regularisation
% algorithm for problems on optimal transport.
%
% **Input:**
% c:   cost matrix of size MxN
% m:   discrete probability vector of size M
% n:   discrete probability vector of size N
% eps: regularisation constant
%
% **Output:**
% x:    best feasible point found after optimisation
% obj:  objective value at x
% temp: time it took to compute x
%
%%
    %Number of iterations
    iters = 1000; 
    
    % Fetch M and N
    M = length(m);
    N = length(n);
    
    % Compute Gibbs Kernel
    Geps = exp(-c/eps);
    
    %Initialisations
    cons  = 1/sum(sum(Geps));
    a     = ones(M,1) * cons;
    b     = ones(N,1) * cons;
    GepsT = Geps.';
    %Initialising the transpose
    
    
    %gam = (b.*((a.*Geps).'));
    %obj = [sum(c*gam,'all')];
    %%initial objective calculation
    
    %starting time.
    tStart = tic;
    for i = 1:iters
        a = m./(Geps*b);
        b = n./(GepsT*a); 
        %updating a and b
        
        %gam = (b.*((a.*Geps).'));
        %obj(end+1) = sum(c*gam,'all');
        %%objective calculation.
        
    end
    % Update time clock
    temp = toc(tStart);
    
    % Final gamma
    gam = (b.*((a.*Geps).')); 
    obj = sum(c .* gam,'all'); 
     
end