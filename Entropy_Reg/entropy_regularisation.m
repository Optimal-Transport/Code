function [gam,obj,temp] = entropy_regularisation(m,n,c,eps)


    %starting time.
    tStart = tic;
    
    iters = 1000; %Number of iterations.
    
    % Fetch lengths of m and n
    N = length(m);
    M = length(n);
    %Save that objective function
    

    
    Geps = exp(-c/eps);
    %Gibbs Kernel
    
    cons = 1/sum(sum(Geps));
    a = ones(M,1)*cons;
    b = ones(N,1)*cons;
    %Initialisations
    
    GepsT = Geps.';
    %Initialising the transpose
    
    
    gam = (b.*((a.*Geps).')).'; 
    obj = [sum(c*gam,'all')];
    %initial objective calculation
    
    
    
    
    

    for i = 1:iters
        a = m./(Geps*b);
        b = n./(GepsT*a); 
        %updating a and b
        
        gam = (b.*((a.*Geps).')).';
        obj(end+1) = sum(c*gam,'all');
        %objective calculation.
        
    end
    
    
    gam = (b.*((a.*Geps).')).';%final gam
    
    
    temp = toc(tStart); %total time
end