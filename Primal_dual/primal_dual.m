function [x, y_1, y_2, obj] = primal_dual(c, m, n, iter)
% Algorithm for calculating solution x, in the primal space
% and y_1, y_2 in the dual space. 
% Also returns the value of the objective function c*x at each
% iteration.

% Initialise tau
tau = 0.5;

% Initialise sigma
sig = 0.5;

% Fetch lengths of m and n
N = length(m);
M = length(n);

% Intialise x
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
