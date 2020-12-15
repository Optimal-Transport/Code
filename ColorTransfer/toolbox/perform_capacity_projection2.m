function x = perform_capacity_projection2(y,mu,somme)

% perform_capacity_projection - projection on the set of valid capacity
% vector according to probability vector mu, such that <x,mu>=somme
%
%    x = perform_capacity_projection(y,mu);
%
%   x is the orthogonal projection of y on the set of x>=0 and
%   <x,mu>=x'*mu=somme
%
%   Copyright (c) 2014 Julien Rabin

Z = size(y); y = y(:); mu = mu(:);
if nargin<2
    error('2 arguments required');
end

n = length(y);

if 0 % old code : bug
    [tmp, idx] = sort(y, 'descend'); % .*mu
    y_sort  = y (idx);
    mu_sort = mu(idx);

    S = cumsum(y_sort .* mu_sort) - somme;
    M = cumsum(mu_sort.^2);
    C = S./M;

    i = 1;
    while ( i<n )
        if (y_sort(i+1) < C(i+1) * mu_sort(i+1))
            break; %error :  we should not stop so early  !
        else
            i=i+1;
        end
    end

    x = max(0, y - C(i)*mu);
    x = reshape(x, Z);
    if abs(mu'*x(:)-somme)>1e-9
        warning('Problem');
    end
    
else
    [tmp, idx] = sort(y./mu, 'descend'); % .*mu
    y_sort  = y (idx);
    mu_sort = mu(idx);

    S = cumsum(y_sort .* mu_sort) - somme;
    M = cumsum(mu_sort.^2);
    C = S./M;

    i = 1;
    while ( i<n )
        if (y_sort(i+1) < C(i+1) * mu_sort(i+1))
            break; %error :  we should not stop so early  !
        else
            i=i+1;
        end
    end

    x = max(0, y - C(i)*mu);
    x = reshape(x, Z);
    if abs(mu'*x(:)-somme)>1e-9
        warning('Problem');
    end
end
