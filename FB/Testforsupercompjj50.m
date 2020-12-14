% Perfomance values
obj  = [];
temp = [];

jj = 50;
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

save('Results/jj50/FB-Temp.mat', 'temp_fb')
save('Results/jj50/FB-Obj.mat', 'obj_fb')


%% Forward-Backwards with backtracking
% Run algorithm
verb = true;  % Return objective evaluation at each point
[x_fbb,obj_fbb,temp_fbb] = FB_Backtracking(c,m,n,verb, 5e-5);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/jj50/FBbacktrack-Temp.mat', 'temp_fbb')
save('Results/jj50/FBbacktrack-Obj.mat', 'obj_fbb')

%% FISTA
% Run algorithm Fista
verb = true;  % Return objective evaluation at each point
%%
[x_fi,obj_fi,temp_fi] = FISTA(c,m,n,verb, 5e-5);

temp(end+1) = temp_fi;
obj(end+1)  = obj_fi(end);

save('Results/jj50/Fista-Temp.mat', 'temp_fi')
save('Results/jj50/Fista-Obj.mat', 'obj_fi')


%% FISTA with backtrack
% Run algorithm
verb = true;  % Return objective evaluation at each point
[x_fib,obj_fib,temp_fib] = FISTA_Backtracking(c,m,n,verb, 5e-5);

temp(end+1) = temp_fib;
obj(end+1)  = obj_fib(end);

save('Results/jj50/Fistawbacktrack-Temp.mat', 'temp_fib')
save('Results/jj50/Fistawbacktrack-Obj.mat', 'obj_fib')

%% Primal-Dual
% Run algorithm
%%
[x_pd, obj_pd, y_1, y_2, temp_pd] = Primal_Dual(c, m, n, true, 5e-5);

temp(end+1) = temp_pd;
obj(end+1)  = obj_pd(end);

save('Results/jj50/Primaldual-Temp.mat', 'temp_pd')
save('Results/jj50/Primaldual-Obj.mat', 'obj_pd')

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

save('Results/jj50/Linearprog-Temp.mat', 'temp_lp')
save('Results/jj50/Linearprog-Obj.mat', 'fval')

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

save('Results/jj50/Entropicreg-Temp.mat', 'temp_er')
save('Results/jj50/Entropicreg-Obj.mat', 'obj_er')

save('Results/jj50/BarChart_temp.mat','temp')
save('Results/jj50/BarChart_finalobj.mat','obj')