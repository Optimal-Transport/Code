function z = prox(y,m,n, tol)
% z = prox(y,m,n,tol) : Computes proximal of the indicator function of the 
% intersection of the simpleces governed by m and n, respectively.
%
% z satisfies sum(y,1) ~ m' and sum(y,2) ~ n with tolerance `tol`.
%
% By default, tol = 1e-10.

    if nargin < 4
        tol = 1e-10;
    end
    
    % Recover size of matrix
    M = length(m);
    N = length(n);
    
    % ? is selected as the midpoint of the interval
    ep = 0.5;
    % The initial points for v_1 and v_2 can be built using the information
    % available already from n and m:
    v_1 = zeros(M,N);
    v_2 = zeros(M,N);
    v_1(:,1) = m;
    v_2(1,:) = n;

    % Define ? as a fixed value
    the = 1.5;
    ith = 1/the;

    % Initialise z:
    z = y - 0.5 * (v_1 + v_2);

    while indicator_simplex(z, m, 2) + indicator_simplex(z, n', 1) ~= 0.0

        % Update v_1 using its projection over simplex C_1^m
        v_1 = v_1 + the * (z - project_simplex(ith * v_1 + z, m, 2) );
        % Update v_2 using its projection over simplex C_2^n
        v_2 = v_2 + the * (z - project_simplex(ith * v_2 + z, n', 1) );
        % Update z
        z = y - 0.5 * (v_1 + v_2);

        % The indicator evaluation is too hard for double precision arithmetic.
        % Thus, we relax this condition:
        if all(abs(sum(abs(z),1) - n') < 1e-10)
            if all(abs(sum(abs(z),2) - m) < 1e-10)
                z = abs(z);
               break 
            end
        end
    end