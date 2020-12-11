function [cost,F] = Optimal_Transport(a,b,C)

% Optimal transport between normalized histograms (i.e. discrete probability density function)
% Warning : this function relies on the 'linprog' Matlab function from
% optimisation toolbox
%
% USAGE :
%           [cost,F] = Optimal_Transport(a,b,C);
%
% SOLVES :
%       min_F <C,F> where C is the [optional] input cost matrix
% subject to the following constraints on the optimal flow matrix F:
%       sum_j F_ij = a_i (sum of mass flows starting from bin i of a is equal to a_i)
%       sum_i F_ij = b_j (                  arriving  at      j    b             b_j)
%       sum_ij F_ij = 1 ( = sum_i a_i = sum_i b_i )
%       1 >= F_ij >= 0 (flows are positives)
%
% INPUTS
%   - a and b are two vectors (normalized histograms) with respectively N
%   and M entries (we do not require N=M)
%   - C is the [optional] N-by-M cost matrix, such that the ground distance
%   between ith and jth bins of histograms a and b respectively is [C_ij].
%   if not specified, the default setting is : C_ij = |i-j|^2 (squared euclidean distance)
%
% OUTPUTS
%   - F is the N-by-M optimal flow matrix where each entry F(i,j) corresponds to the mass transportation from i to j of the amount of a(i) into b(j)
%   - cost = <C,F> is the total transportation cost (so called L2-Wasserstein squared distance with default settings)
%
% Julien.Rabin@unicaen.fr 2013 (c) GREYC, Université de Caen
%
% % [example 1] with L1 ground cost
% N = 1e2; x = linspace(0,1,N);
% gauss = @(x,m,s) exp(-.5*(x-m).^2/s^2 );
% a = gauss(x,.3,.05); b = gauss(x,.6,.15);
% a = a/sum(a); b = b/sum(b);
% figure, bar(a,'b'), hold on, bar(-b,'r')
% t = (0:1:N-1)/N; C = toeplitz(t); % L1 distance matrix
% [cost,F] = Optimal_Transport(a,b,C);
% % figure, imagesc(F), title('optimal flow from a to b')
% W2_L1 = sum(abs(cumsum(a-b)))/N % for sanity check (special case with known as L1-Wasserstein distance)
% cost, % should be equal to W2_L1
%
% % [example 1 bis] same example with L2^2 ground cost 
% C = toeplitz(t.^2); % L2^2 distance matrix
% [cost,F] = Optimal_Transport(a,b,C);
% figure, imagesc(F), title('optimal flow from a to b')
%
% % [example 2]
% N = 5; M = 10;
% a = 1:N; b = M:-1:1;
% a = a/sum(a); b = b/sum(b);
% figure, bar(a,'b'), hold on, bar(-b,'r')
% x = linspace(0,1,N); y = linspace(0,1,M); C = zeros(N,M);
% for i=1:N, for j=1:M, C(i,j) = (x(i)-y(j)).^2; end, end % L2 distance matrix
% [cost,F] = Optimal_Transport(a,b,C);
% figure, imagesc(F), title('optimal flow from a to b')
%
% % [example 3] 2D optimal transport
% n = 10; N = n*n; x = linspace(0,1,n);
% [X,Y] = meshgrid(x,x);
% gauss = @(x,m,s) exp(-.5*(((x-m)./s).^2));
% a = gauss(X,.3,.05).*gauss(Y,.2,.1); b = gauss(X,.6,.15).*gauss(Y,.7,.07);
% a = a./sum(a(:)); b = b./sum(b(:)); figure, imagesc(a-b), title('mass distribution')
% x = [X(:), Y(:)]; y = repmat( sum(x.^2,2), [1 N]);
% C = y + y' - 2*x*x'; figure, imagesc(C), title('distance cost matrix')
% [cost,F] = Optimal_Transport(a(:),b(:),C);
% figure, imagesc([a,inf(n,1),b]), title('optimal flow from a to b')
% K = 10; [val,idx] = sort(F(:),'descend');
% idx = idx(1:K);
% for I=idx(:)'
%     [i,j,k,l] = ind2sub([n n n n],I);
%     h = line([j;l+n+1],[i,k],'color','g');
% end
%
% See Also: Optimal_Transport, Optimal_Relaxed_Transport

%% 1 - check inputs
if nargin<2
    disp('error: not enough inputs !')
    return
% elseif length(a)~=length(b)
%     disp('error: input histograms must have the same size !')
%     return
elseif max(size(C)~=[length(a) length(b)])
    disp('error: cost matrix must have the same size than input histograms!')
    return
end

N = length(a); a = a(:);
M = length(b); b = b(:);
NM = N*M;

if nargin<3
    if N~=M
        disp('error: input histograms must have the same size when no cost matrix is given!')
        return
    end
    disp('Default: cost matrix is the squared euclidean distance matrix !')
    t = (0:1:N-1); C = toeplitz(t.^2);
end

if abs(sum(a)-1)>1e-15
  disp('Warning: normalization of a !')
  a = a/sum(a);
end
if abs(sum(b)-1)>1e-15
  disp('Warning: normalization of b !')
  b = b/sum(b);
end

c = reshape(C, NM, 1); %lecture en colonne

%% 2 - build constraint matrices for flow f, viewed as [N*M by 1] array (read column-wise)
disp('Build constraint matrix')
tic
% a) upper and lower bound on F: ub >= f >= lb
% each flow is positive and less than 1
lb = zeros(NM,1); 
ub = []; % ones(NM,1);

% b) affine inequalities: Af < B
A=[]; B=[]; % no inequalities

% c) affine equalities: Aeq f = beq
% using SPARSE representation
beq = cat(1,a,b); % [N+M by 1] array

% A1 * f = a : equivalent to F * ones(M,1) = a
A1 = repmat(speye(N),1,M);
% A2 * f = b : equivalent to ones(1,N) * F = b
A2 = sparse(M,NM);
for i=1:M
    A2(i,(i-1)*N+1:i*N)=ones(1,N);
end

% adding sum(f) = 1 constraint is unnecessary

Aeq = cat(1,A1,A2);
clear('A1','A2');
toc


%% 3 - optimisation
%opt = optimset;
opt = optimset('Display','off','Maxiter',1000);
tic
    [f,cost,EXITFLAG] = linprog(c,A,B,Aeq,beq,lb,ub,[],opt);
toc
if EXITFLAG~=1
    display(['Problem with linprog. Exitflag=' num2str(EXITFLAG)]);
end
%% 4 - set outputs

F = reshape(f,N,M);

% verif
% F *ones(M,1) - a
% F'*ones(N,1) - b

