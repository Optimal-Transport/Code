function [x,obj,temp] = Forward_Backward(c,m,n,tol)
%%
% Forward_Backward(c,m,n,tol) : Executes the Forward-Backward algorithm
% applied to problems on optimal transport.
%
% **Input:**
% c:   cost matrix of size MxN
% m:   discrete probability vector of size M
% n:   discrete probability vector of size N
% tol: numerical tolerance of the algorithm: it stops if the norm between a
%      pair of iterations is less than this value (default tol = 1e-4)
%      
% **Output:**
% x:    best feasible point found after optimisation
% obj:  objective value at x
% temp: time it took to compute x
%%

    if nargin < 6
        tol = 1e-4;
    end
    % Recover M and N
    M = length(m);
    N = length(n);
    
    % First compute ?
    mu = 1000;%norm(c,2);        %% 1 -> 10^-1 -> 10^-2 -> ...
    % ? is selected as the midpoint of the interval
    e = 1/mu;%0.5 * 1/mu;        %% remove
    % ? does not depend on the current iteration
    gam = 2/mu - e;            %% 1.9/mu
    % likewise, we do not require a change in ?
    lam = 0.5 * (1 + e);   %% 1.0

    %% 
    % x_0 is projected to be a feasible initial point 
    %%
    %x = (m' + n)/(N+M);
    x = zeros(M,N);
    x(:,1) = m;
    x(1,:) = n;
    x(1,1) = 0.5 * (m(end) + n(end));
    % The initial points for v_1 and v_2 can be built using the information
    % available already from n and m:
    v_1 = zeros(M,N);
    v_2 = zeros(M,N);
    v_1(:,1) = m;
    v_2(1,:) = n;
    % Compute proximal operator at C and update v_1 and v_2
    [x, v_1, v_2] = prox_i(x,m,n,v_1,v_2);
    
    % iters controls the number of iterations
    iters = 100;                                 %% increase
    % The distance between points will serve as stopping criteria
    norm_difference = Inf;
    % Objective value
    obj = Inf;
    % Measure time
    tStart = tic;
    
    %% Now we perform the FB iteration:
    for it = 1:iters
        y = x - gam * c;                            % gam = 10 takes > 10 min for second instance
        % Proximal operation
        [y, v_1, v_2] = prox_i(y,m,n,v_1,v_2);
        u = (1-lam) * x + lam * y;
        norm_difference = norm(x-u);
        x = u;
        if norm_difference < tol * norm(u)
            break
        end
    end
    % Update time clock
    temp = toc(tStart);
    obj  = sum(sum(c.*x));
    
    % See order of magnitude and number of iterations
    [log(norm_difference)/log(10), it]  
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    