jj=4;

obj  = [];
temp = [];
M = 200 + jj*200;          N = 200 + jj*200;
m = (1/M) * ones(M,1);     n = (1/N) * ones(N,1);

X = linspace(0,1,M)';
Y = linspace(1,2,M)';

c = squareform(pdist([X;Y],'squaredeuclidean'));
c = c(1:M, M+1:end);

clear X Y;

%Lin progr system Ax = b
A = [kron( ones(1,N), speye(M) ); kron( speye(N), ones(1,M) )];
b = [ m;n ];

    %Lin progr system Ax = b
    A = [kron( ones(1,N), speye(M) ); kron( speye(N), ones(1,M) )];
    b = [ m;n ];

%Run algorithm
options = optimoptions('linprog','Algorithm','dual-simplex','OptimalityTolerance', 0.1, 'ConstraintTolerance',1e-4);
tStart = tic;
[x_lp,fval,exitflag,output,lambda] = linprog(c(:),[],[],A,b, sparse(zeros(N*M,1)), [], options);

temp_lp = toc(tStart);
temp(end+1) = temp_lp;
fval_sum = sum(fval);
obj(end+1) = sum(fval);

save('Results/Linprog/Linearprog-Temp5e-01.mat', 'temp_lp')
save('Results/Linprog/LinprogLinearprog-fval5e-01.mat', 'fval_sum')

for kk = [2:5]
    %Lin progr system Ax = b
    A = [kron( ones(1,N), speye(M) ); kron( speye(N), ones(1,M) )];
    b = [ m;n ];

    %Run algorithm
    options = optimoptions('linprog','Algorithm','dual-simplex','OptimalityTolerance', 5e-kk, 'ConstraintTolerance',1e-6);
    tStart = tic;
    [x_lp,fval,exitflag,output,lambda] = linprog(c(:),[],[],A,b, sparse(zeros(N*M,1)), [], options);

    temp_lp = toc(tStart);
    temp(end+1) = temp_lp;
    fval_sum = sum(fval);
    obj(end+1) = sum(fval);

    save('Results/Linprog/Linearprog-Temp5e-' + string(kk) + '.mat', 'temp_lp')
    save('Results/Linprog/LinprogLinearprog-fval5e-' + string(kk) + '.mat', 'fval_sum')
end
save('Results/Linprog/Linprog_Tempcrit.mat', 'temp')
