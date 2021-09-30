jj =4;
% Perfomance values
obj  = [];
temp = [];
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
[x_fb,obj_fb,temp_fb, tempcrit_fb] = Forward_Backward(c,m,n,verb, 5e-5);

temp(end+1) = temp_fb;
obj(end+1)  = obj_fb(end);

save('Results/jj' + string(jj) + '/FB-Tempcrit_jj=' + string(jj) + '.mat', 'tempcrit_fb')
save('Results/jj' + string(jj) + '/FB-Tempjj=' + string(jj) + '.mat', 'temp_fb')
save('Results/jj' + string(jj) + '/FB-Objjj=' + string(jj) + '.mat', 'obj_fb')


%% Forward-Backwards with backtracking
% Run algorithm
verb = true;  % Return objective evaluation at each point
[x_fbb,obj_fbb,temp_fbb, tempcrit_fbb] = FB_Backtracking(c,m,n,verb, 5e-5);

temp(end+1) = temp_fbb;
obj(end+1)  = obj_fbb(end);

save('Results/jj' + string(jj) + '/FBbacktrack-Tempcritjj' + string(jj) + '.mat', 'tempcrit_fbb')
save('Results/jj' + string(jj) + '/FBbacktrack-Temp.matjj' + string(jj) + '.mat', 'temp_fbb')
save('Results/jj' + string(jj) + '/FBbacktrack-Obj.matjj' + string(jj) + '.mat', 'obj_fbb')

%% FISTA
% Run algorithm Fista
verb = true;  % Return objective evaluation at each point
%%
[x_fi,obj_fi,temp_fi, tempcrit_fi] = FISTA(c,m,n,verb, 5e-5);

temp(end+1) = temp_fi;
obj(end+1)  = obj_fi(end);

save('Results/jj' + string(jj) + '/Fista-Tempcritjj=' + string(jj) + '.mat', 'tempcrit_fi')
save('Results/jj' + string(jj) + '/Fista-Tempjj=' + string(jj) + '.mat', 'temp_fi')
save('Results/jj' + string(jj) + '/Fista-Objjj=' + string(jj) + '.mat', 'obj_fi')


%% FISTA with backtrack
% Run algorithm
verb = true;  % Return objective evaluation at each point
[x_fib,obj_fib,temp_fib,tempcrit_fib] = FISTA_Backtracking(c,m,n,verb, 5e-5);

temp(end+1) = temp_fib;
obj(end+1)  = obj_fib(end);

save('Results/jj' + string(jj) + '/Fistawbacktrack-Tempcritjj=' + string(jj) + '.mat', 'tempcrit_fib')
save('Results/jj' + string(jj) + '/Fistawbacktrack-Tempjj=' + string(jj) + '.mat', 'temp_fib')
save('Results/jj' + string(jj) + '/Fistawbacktrack-Objjj=' + string(jj) + '.mat', 'obj_fib')

%% Primal-Dual
% Run algorithm
%%
[x_pd, obj_pd, y_1, y_2, temp_pd,tempcrit_pd] = Primal_Dual(c, m, n, true, 5e-5);

temp(end+1) = temp_pd;
obj(end+1)  = obj_pd(end);

save('Results/jj' + string(jj) + '/Primaldual-Tempcritjj=' + string(jj) + '.mat', 'tempcrit_pd')
save('Results/jj' + string(jj) + '/Primaldual-Tempjj=' + string(jj) + '.mat', 'temp_pd')
save('Results/jj' + string(jj) + '/Primaldual-Objjj=' + string(jj) + '.mat', 'obj_pd')

%%%%%%%%%%%%%%%
% Lin progr system Ax = b
%A = [kron( ones(1,N), speye(M) ); kron( speye(N), ones(1,M) )];
%b = [ m;n ];

% Run algorithm
%options = optimoptions('linprog','Algorithm','dual-simplex','OptimalityTolerance', 5e-5, 'ConstraintTolerance',1e-6);
%tStart = tic;
%[x_lp,fval,exitflag,output,lambda] = linprog(c(:),[],[],A,b, sparse(zeros(N*M,1)), [], options);

%temp_lp = toc(tStart);
%temp(end+1) = temp_lp;
%obj(end+1) = sum(fval);

%save('Results/Linearprog-Temp.mat', 'temp_lp')
%save('Results/Linearprog-fval.mat', 'obj(end)')

%%%%%%%%%%%%%%%%%%%
% Perfomance values Entropic- Reg
eps = 0.0025; % Regularisation constant
verb = true;  % Return objective evaluation at each point
%% 
% Run algorithm
%%
[x_er,obj_er,temp_er, tempcrit_er] = Entropy_Regularisation(m,n,c,eps,verb, 5e-5);
temp(end+1) = temp_er;
obj(end+1)  = obj_er(end);

save('Results/jj' + string(jj) + '/Entropicreg-Tempcritjj=' + string(jj) + '.mat', 'tempcrit_er')
save('Results/jj' + string(jj) + '/Entropicreg-Tempjj=' + string(jj) + '.mat', 'temp_er')
save('Results/jj' + string(jj) + '/Entropicreg-Objjj=' + string(jj) + '.mat', 'obj_er')