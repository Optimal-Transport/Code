% Perfomance values
obj  = [];
temp = [];
eps = 0.01;

for jj = 0:5
    M = 200 + jj*200;                   N = 100 + jj*100;
    m = (1/M) * ones(M,1);     n = (1/N) * ones(N,1);
    
    for ii = 1:2
        % Compute distances
        rng(2020)
        X = rand(M,1);
        Y = 1 + rand(N,1);
        if      ii == 1
            c = squareform(pdist([X;Y]));
            c = c(1:M, M+1:end);
            
        elseif  ii == 2
            c = squareform(pdist([X;Y],'squaredeuclidean'));
            c = c(1:M, M+1:end);
            
        else
            c = squareform(pdist([X;Y],'minkowski',0.5));
            c = c(1:M, M+1:end);
            
        end
        clear X Y;
        
        %% 
        % Run algorithm
        %%
        [gam,objm,tempm] = entropy_regularisation(m,n,c,eps);
        temp(end+1) = tempm;
        obj(end+1) = objm;
    end
end