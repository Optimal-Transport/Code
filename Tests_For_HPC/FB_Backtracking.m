function [x,obj,temp,temp_crit] = FB_Backtracking(c,m,n,collect_obj,tol)
%%
% FB_Backtracking(c,m,n,collect,tol) : Executes the Forward-Backward algorithm
% applied to problems on optimal transport with backtracking.
%
% **Input:**
% c:   cost matrix of size MxN
% m:   discrete probability vector of size M
% n:   discrete probability vector of size N
% collect_obj: boolean value; if true, then all objective values are stored
% tol: numerical tolerance of the algorithm: it stops if the norm between a
%      pair of iterations is less than this value (default tol = 1e-4)
%
% **Output:**
% x:    best feasible point found after optimisation
% obj:  objective value at x
% temp: time it took to compute x
%%

    if nargin < 5
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
    iters = 1000;                                 %% increase
    % The distance between points will serve as stopping criteria
    norm_difference = Inf;
    % Objective value
    obj = [];
    % initial objective calculation
    if collect_obj
        obj = [sum(c.*x,'all')];
    end
    
    % Measure time
    tStart = tic;
    
    % Display time at different intervals
    temp_crit = [];
    crit_tol = [];
    for i = [0:5]
        crit_tol(end+1) = 5 * 10 ^ (-i);
        crit_tol(end+1) = 1 * 10 ^ (-i);
    end
    j = 1;

    %% Now we perform the FB iteration:
    for it = 1:iters
        y = x - gam * c;                            % gam = 10 takes > 10 min for second instance
        % Proximal operation
        [y, v_1, v_2] = prox_i(y,m,n,v_1,v_2);
        
        while 4 * sum(c.*(y-x), 'all') > mu * sumsqr( y - x )
            mu  = mu * 2;
            gam = 0.5 * gam;
            y = x - gam * c;
            [y, v_1, v_2] = prox_i(y,m,n,v_1,v_2);
        end
        
        u = (1-lam) * x + lam * y;
        norm_difference = norm(x-u);
        
        %Average objective from 2 previous
        aver_obj = (sum(c.*u,'all') + (sum(c.*x,'all')))/2;
        
        x = u;
        
        % Store objective if needed
        if collect_obj
            obj(end+1) = sum(sum(c.*x));
        end
        
        
        % Store temp for certain tolerance
        if abs(aver_obj - 1) < crit_tol(j)
            temp_crit(end+1) = toc(tStart);
            j = j + 1;
        end
        % Check tolerance
        if abs(aver_obj - 1) < tol
            break
        end
    end
    % Update time clock
    temp = toc(tStart);
    temp_crit(end+1) = temp;
    obj(end+1) = sum(sum(c.*x));

    % See order of magnitude and number of iterations
    [log(norm_difference)/log(10), it]
