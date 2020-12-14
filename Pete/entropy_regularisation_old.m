function [x,obj,temp] = entropy_regularisation_old(m,n,c,eps,collect_obj,tol)


    %starting time.
    tStart = tic;
    
    iters = 10000; %Number of iterations.
    
    % Fetch lengths of m and n
    N = length(n);
    M = length(m);
    %Save that objective function
    
    obj = [];
    
    Geps = exp(-c/eps);
    %Gibbs Kernel
    %disp(size(Geps))
    
    cons = 1/sum(sum(Geps));
    a = ones(M,1)*cons;
    b = ones(N,1)*cons;
    %disp(size(a))
    %disp(size(b))
    %Initialisations
    
    GepsT = Geps.';
    %Initialising the transpose
    
    
    %gam = (b.*((a.*Geps).'));
    %obj = [sum(c*gam,'all')];
    %%initial objective calculation
    
    
    
    
    

    for i = 1:iters
        a = m./(Geps*b);
        b = n./(GepsT*a); 
        %updating a and b
        
        %gam = (b.*((a.*Geps).'));
        %obj(end+1) = sum(c*gam,'all');
        %%objective calculation.
        %Calculate x
        x = (b.*((a.*Geps).'));
        % Store objective if needed
        if collect_obj
            obj(end+1) = sum(sum(c.*x));
        end
        if abs(sum(sum(c.*x)) - 1) < tol
            break
        end
        
    end
    %calculating final objective
    
   
    
    
    temp = toc(tStart); %total time
end