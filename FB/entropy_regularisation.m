function [x,obj,temp] = Entropy_Regularisation(m,n,c,eps,collect_obj)
%%
% Entropy_Regularisation(m,n,c,eps) : Executes the Entropy Regularisation
% algorithm for problems on optimal transport.
%
% **Input:**
% c:   cost matrix of size MxN
% m:   discrete probability vector of size M
% n:   discrete probability vector of size N
% eps: regularisation constant (so far the best value has been 0.0025)
% collect_obj: boolean value; if true, then all objective values are stored
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
    
    % Initialisations
    cons  = 1/sum(sum(Geps));
    a     = ones(M,1) * cons;
    b     = ones(N,1) * cons;
    
    obj = [];
    % initial objective calculation
    if collect_obj
        gam = a.*Geps .* b';
        obj = [sum(c.*gam,'all')];
    end
    
    % Starting time
    tStart = tic;
    for i = 1:iters
        % Updating a and b
        a = m./(Geps*b);
        b = n./(a'*Geps)';
        
        if collect_obj
            gam = a.*Geps .* b';
            obj(end+1) = sum(c.*gam,'all');
        end
        
    end
    % Update time clock
    temp = toc(tStart);
    % Final estimation
    x   = a.*Geps .* b';
    obj(end+1) = sum(c .* x,'all'); 
    
end