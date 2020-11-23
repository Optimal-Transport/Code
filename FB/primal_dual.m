function [x, obj, y_1, y_2, temp] = Primal_Dual(c, m, n, tol)
%%
% Primal_Dual(c, m, n, tol) : Executes the Primal-Dual algorithm
% applied to problems on optimal transport.
%
% **Input:**
% c:   cost matrix of size MxN
% m:   discrete probability vector of size M
% n:   discrete probability vector of size N
% tol: numerical tolerance of the algorithm: it stops if the norm between a
%      pair of iterations is less than this value (default tol = 1e-10)
%      
% **Output:**
% x:    best feasible point found after optimisation
% obj:  objective value at x
% y:    Pair of dual variables
% temp: time it took to compute x
%%

    if nargin < 6
        tol = 1e-4;
    end
    % Recover M and N
    M = length(m);
    N = length(n);

    % Initialise tau
    tau = 0.5;
    % Initialise sigma
    sig = 0.5;

    %% x_0 is projected to be a feasible initial point 
    x = (m' + n)/(N + M);
    
    
% Intialise y_1 and y_2
y_1 = zeros(N,M);
y_2 = zeros(N,M);
y_1(1,:) = m;
y_2(:,1) = n;

%Save that objective function
obj = [sum(c*x,'all')]

for k = 1:iter
    % Update x
    xnew = x - tau * c - tau *(y_1 + y_2);
    % Update x using the projection over the simplex C_1^m
    y_1new = y_1 + sig * (2 * xnew - x) - sig * project_simplex((1/sig) * y_1 + 2 * xnew - x, m, 2);
    % Update x using the projection over the simplex C_2^n
    y_2new = y_2 + sig * (2 * xnew - x) - sig * project_simplex((1/sig) * y_2 + 2 * xnew - x, n', 1);
        
    % Reset x, y_1, y_2 for the next iteration
    x = xnew;
    y_1 = y_1new;
    y_2 = y_2new;

    % Update objective function
    obj(end+1) = sum(c*x,'all');
    
end
end
