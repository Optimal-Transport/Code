function [D,P] = sinkhorn_distance(a,b,C,lambda,options)
% compute sinkhorn distance based on 
% 'lightspeed optimal transport' by Marco Cuturi 2013
%
% inputs :  - a [N x 1] & b  [M x 1] are two normalized histograms such that sum(a) = sum(b) = 1
%           - C [N x M] is a cost matrix of affectation of bins of a to bins of j
%           - lambda is a regularization parameter, so that : 
%               * lambda = 0 leads to maximum entropy transport (every bin is transported to every bin)
%               % lambda = infty solves the optimal transport problem
%           - options is a structure with miscellanous options such as
%           niter (number of iterations)
%
% outputs : - D is the sinkhorn distance : <P,C>
%           - P is the optimal solution of : min_{P \in S} < P, C + log(P) /lambda >
%             where S is the set of bistochastic matrices : P>0, P1=a, P'1=b
%
% exemple : 
%     clear all, close all, clc
%     N = 400; M = 500;
%     gauss = @(x,m,s) exp(-.5*(x-m).^2/s^2 )/s/sqrt(2*pi);
%     x = linspace(0,1,N);
%     y = linspace(0,1,M);
%     t = 1; y = y+t; % Y is translated
%     % 2 gaussians
%     a = .8*gauss(x,.15,.05)  + .2*gauss(x, .55,.05);  
%     b = .9*gauss(y,t+.9,.05) + .5*gauss(y,t+.1,.02);
%     a = a(:) / sum(a); b = b(:) / sum(b);
%     X = repmat(x(:), [1,M]); Y = repmat(y(:)', [N,1]);
%     C = X.^2 + Y.^2 - 2 * X.*Y; % figure, imagesc(C) % quadratic (L_2^2) euclidean cost matrix
%     C_L2 = C / max(C(:)); % normalization in [0,1]
%
%     options.niter = 1e2;
%     [D,P] = sinkhorn_distance(a,b,C_L2,1e2,options);
% 
%     figure, plot(x,a,'r'), hold on, plot(y,b,'b')
%     figure, 
%     subplot(1,2,1), imagesc(a(:)*b(:)'), title('lambda = \infty')
%     subplot(1,2,2), imagesc(P), title(['lambda = ' num2str(lambda)])
%
% J. Rabin (c) 2014 GREYC julien.rabin@unicaen.fr

if lambda < 0
    disp('error : regularization paramter should be positive ')
    return
end

if exist('options')
    niter = options.niter;
else
    disp('default number of iterations : 100')
    niter = 1e2;
end

a = a(:);
b = b(:);

N = numel(a);
M = numel(b);
if size(C,1) ~= N || size(C,2) ~= M
    disp('error : cost matrix does not fit input vectors')
    return
end

K = exp(-lambda*C); % point-wise exponential
Ka = bsxfun(@rdivide, K , a); % equivalent to diag(1./a)*K for row normalization
Kb = bsxfun(@rdivide, K', b);


tic
u = ones(N,1)/N; % init
for it = 1 : niter % algorithm
    u = 1./( Ka * ( 1 ./ (Kb * u) ) );
end
toc

v = 1 ./ (Kb * u);
P = bsxfun(@times, v, (bsxfun(@times, K, u))')'; % equivalent to : diag(u) * K * diag(v); % 
% optimal and regularized matrix
D = sum( P(:) .* C(:) ); % sinkhorn transport cost distance

