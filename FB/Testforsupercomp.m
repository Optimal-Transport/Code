% Perfomance values
obj  = [];
temp = [];

jj = 5;
M = 200 + jj*200;          N = 200 + jj*200;
m = (1/M) * ones(M,1);     n = (1/N) * ones(N,1);

X = linspace(0,1,M)';
Y = linspace(1,2,M)';
        
c = squareform(pdist([X;Y],'squaredeuclidean'));
c = c(1:M, M+1:end);
%save('Instances/Cost_SE_' + string(M) + '_' + string(N) + '.mat', 'c')  % Store matrix for reproducibility
clear X Y;

%% Forwards backwards
% Run algorithm
verb = true;  % Return objective evaluation at each point
%%
[x_fb,obj_fb,temp_fb] = Forward_Backward(c,m,n,verb, 5e-5);

temp(end+1) = temp_fb;
obj(end+1)  = obj_fb(end);

save('Results/FB-Temp.mat', 'temp')
save('Results/FB-Obj.mat', 'obj')


%% Forward-Backwards with backtracking
% Run algorithm
verb = true;  % Return objective evaluation at each point
[x_fbb,obj_fbb,temp_fbb] = FB_Backtracking(c,m,n,verb, 5e-5);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/FBbacktrack-Temp.mat', 'temp')
save('Results/FBbacktrack-Obj.mat', 'obj')

%% FISTA
% Run algorithm Fista
verb = true;  % Return objective evaluation at each point
%%
[x_fi,obj_fi,temp_fi] = FISTA(c,m,n,verb, 5e-5);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/Fista-Temp.mat', 'temp')
save('Results/Fista-Obj.mat', 'obj')


%% FISTA with backtrack
% Run algorithm
verb = true;  % Return objective evaluation at each point
[x_fib,obj_fib,temp_fib] = FISTA_Backtracking(c,m,n,verb, 5e-5);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/Fistawbacktrack-Temp.mat', 'temp')
save('Results/Fistawbacktrack-Obj.mat', 'obj')

%% Primal-Dual
% Run algorithm
%%
[x_pd, obj_pd, y_1, y_2, temp_pd] = Primal_Dual(c, m, n, true, 5e-5);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/Primaldual-Temp.mat', 'temp')
save('Results/Primaldual-Obj.mat', 'obj')

%%%%%%%%%%%%%%%
% Lin progr system Ax = b
A = [kron( ones(1,N), speye(M) ); kron( speye(N), ones(1,M) )];
b = [ m;n ];
    
% Run algorithm
options = optimoptions('linprog','Algorithm','dual-simplex','OptimalityTolerance', 5e-5, 'ConstraintTolerance',1e-6);
tStart = tic;
[x_lp,fval,exitflag,output,lambda] = linprog(c(:),[],[],A,b, sparse(zeros(N*M,1)), [], options);

temp_lp = toc(tStart);
temp(end+1) = temp_lp;
obj(end+1) = sum(fval);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/Linearprog-Temp.mat', 'temp')
save('Results/Linearprog-Obj.mat', 'obj')

%%%%%%%%%%%%%%%%%%%
% Perfomance values Entropic- Reg
eps = 0.0025; % Regularisation constant
verb = true;  % Return objective evaluation at each point
%% 
% Run algorithm
%%
[x_er,obj_er,temp_er] = Entropy_Regularisation(m,n,c,eps,verb);
temp(end+1) = temp_er;
obj(end+1)  = obj_er(end);

save('Results/Entropicreg-Temp.mat', 'temp')
save('Results/Entropicreg-Obj.mat', 'obj')
