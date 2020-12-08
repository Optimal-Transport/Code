
% This code is free to use for any non-commercial purposes.
%
% Color transfer between RGB images
%
%If using please cite the 2 papers:

% Non-convex relaxation of optimal transport for color transfer, NIPS
% Workshop 2014, Rabin and Papadakis
%
% and
%
% SCALP: Superpixels with Contour Adherence using Linear Path, ICPR 2016,
% Giraud, Ta, Papadakis
%
% (C) Nicolas Papadakis and Julien Rabin, 2017

clear all
close all

% require linprog to compute true Optimal transport for initialization
% if not available, initialization with sinkhorn algorithm
USE_SINKHORN=0;
testlinprog=exist('linprog');  %better with linprog...
if testlinprog ~=2
    disp('WARNING: linprog not available: intialization with sinkhorn giving poor results');
    USE_SINKHORN=1;
end

addpath('./toolbox/')
rep = './images/'

% image selection

name1 = 'fleur-1.jpg', name2 = 'fleur-2.jpg'
name1 = 'flower_a.jpg', name2 = 'flower_b.jpg'

%name1 = 'parrot-1.jpg', name2 = 'parrot-2.jpg'
name1 = 'im5.jpg', name2 = 'im6.jpg'
name1 = 'paris.png', name2 = 'plage2.jpg'
% load image
[x_,x_map] = imread([rep, name1]);
grey=0;
if size(x_,3)==1
    grey=1;
    x_(:,:,2)=x_(:,:,1);
    x_(:,:,3)=x_(:,:,1);
end
% load image

[y_,x_map] = imread([rep, name2]);

% % resize image
% tx = max(size(x_));
% ty = max(size(y_));
% %t = max(cat(2,size(x)),size(y));
% t_max = 2^12; % 2^9
% if tx > t_max
%     k = t_max/tx;
%     x_ = imresize(x_,k);
%     % y = imresize(y,k);
% end
% if ty > t_max
%     k = t_max/ty;
%     y_ = imresize(y_,k);
% end

%% Parameter setting

SHOW=1;
%choose method


KLUSTER = 121;
USE_FAST_SYNTHESIS=1;


ALPHA = 1e3; % 1e2; for cubic formulation of variance, 1e3 for quadratic forumlation
LAMBDA= 4e2; % 4e2: for regularization
RHO = 1e+2; % 1e2: for fidelity term

NB_NEIGHBORS = 30;
VARIANCE_COLOR=1e-1;
VARIANCE_SPATIAL=1e1;
do_yreg=0;
ITMAX=1000;

if 1
    boxname='Set options and parameters ';
    prompt={'Nb of clusters: (integer expected >=1)',...
        '{\bf Color variance parameter} :',...
        'Regularization parameter (>=0):', 'Fidelity parameter of color palette',...
        'Number of neighbors for regularization (>0):', ... % 'Minimum mass to reach on target''s clusters (<=1):',     'Maximum mass to reach on target''s clusters (>=1):', ...
        'Color variance for cluster weights','Spatial variance for cluster weights',...
        'Color regularization for Y: (0: No, otherwise: Yes)',...
        'Use fast synthesis: (set 1 if you have a lot of RAM)',...
        'Number ot iterations : '
        };
    numlines=1;
    defaultanswer={num2str(KLUSTER),  ...
        num2str(ALPHA),num2str(LAMBDA),num2str(RHO),...
        num2str(NB_NEIGHBORS), ...
        num2str(VARIANCE_COLOR), num2str(VARIANCE_SPATIAL), num2str(do_yreg), num2str(USE_FAST_SYNTHESIS), num2str(ITMAX)};
    
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    
    answer=inputdlg(prompt,boxname,numlines,defaultanswer,options);
    
    i = 1;
    KLUSTER         = max(1,str2num(answer{i})); i = i+1;
    KLUSTER =double(uint16(sqrt(KLUSTER)))^2;
    ALPHA = str2num(answer{i}); i = i+1;
    LAMBDA = max(0,str2num(answer{i})); i = i+1;
    RHO = str2num(answer{i}); i = i+1;
    NB_NEIGHBORS = max(1,str2num(answer{i})); i = i+1;
    VARIANCE_COLOR = max(0,str2num(answer{i})); i = i+1;
    VARIANCE_SPATIAL = max(0,str2num(answer{i})); i = i+1;
    
    do_yreg = max(0,str2num(answer{i})); i = i+1;
    
    USE_FAST_SYNTHESIS=str2num(answer{i}); i = i+1;
    ITMAX = max(0,str2num(answer{i})); i = i+1;
    
    
    
    
    
    
end
drawnow;

%% image pre-processing

%%%%%%%%%%%%%%%%%%
%1) sub-sampling %
%%%%%%%%%%%%%%%%%%
x0=x_;
y0=y_;


Nx_X=size(x_,1);
Ny_X=size(x_,2);

Nx_Y=size(y_,1);
Ny_Y=size(y_,2);


spatial_weight = 255;
[YYX,XXX] = meshgrid(linspace(0,spatial_weight,Ny_X), linspace(0,spatial_weight,Nx_X)); % spatial coordonates are considered as colors in [0;255]
XXX = reshape(XXX, [1 Nx_X*Ny_X]);
YYX = reshape(YYX, [1 Nx_X*Ny_X]);


[YYY,XXY] = meshgrid(linspace(0,spatial_weight,Ny_Y), linspace(0,spatial_weight,Nx_Y));
XXY = reshape(XXY, [1 Nx_Y*Ny_Y]);
YYY = reshape(YYY, [1 Nx_Y*Ny_Y]);


if(SHOW)
    figure; subplot(1,2,1); imshow(uint8(x_),[0 1]); title(['source image: [' num2str(size(x_)) ']']);%color images are floats in [0,1]
    subplot(1,2,2); imshow(uint8(y_),[0 1]); title(['target image: [' num2str(size(y_)) ']']);drawnow
end





%Reshape
X_ =zeros(3,Nx_X*Ny_X);
X_(1,:) = reshape(x_(:,:,1) ,1,Nx_X*Ny_X);
X_(2,:) = reshape(x_(:,:,2) ,1,Nx_X*Ny_X);
X_(3,:) = reshape(x_(:,:,3) ,1,Nx_X*Ny_X);
X_= double(X_);

Y_ =zeros(3,Nx_Y*Ny_Y);
Y_(1,:) = reshape(y_(:,:,1) ,1,Nx_Y*Ny_Y);
Y_(2,:) = reshape(y_(:,:,2) ,1,Nx_Y*Ny_Y);
Y_(3,:) = reshape(y_(:,:,3) ,1,Nx_Y*Ny_Y);
Y_= double(Y_);



X_ = [X_; XXX; YYX];
Y_ = [Y_; XXY; YYY];


%% Let's apply the color transfer algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 : perform super-pixel clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx = length(X_);
Ny = length(Y_);

K = KLUSTER; %nombre de clusters


tic
if K>1
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % SUPERPIXELS
    %%%%%%%%%%%%%%%%%%%%%%%%
    %
    
    
    [Kx,Mu_x] = SCALP_segmentation(single(x_),K); disp('1 of 2')
    
    [Ky,Mu_y] = SCALP_segmentation(single(y_),K); disp('2 of 2')
    KKX=max(Kx(:));
    KKY=max(Ky(:));
    
    
    
    Mu_x=Mu_x(1:KKX,:);
    Mu_y=Mu_y(1:KKY,:);
    
    
    Mu_x=Mu_x'; Mu_y=Mu_y'; % (3 or 5)xK %matrix with K clusters
    Kx=Kx'; Ky=Ky'; %1xN matrix of assignation of K clusters 'C'
    
elseif K==1
    Mu_x = sum(X_,2)./Nx;
    Mu_y = sum(Y_,2)./Ny;
    Kx = ones(1,Nx);
    Ky = ones(1,Ny);
else
    return
end

toc;

disp('Creation of cluster information');


P_x = zeros(1,KKX);
P_y = zeros(1,KKY);

P_x = sum( (ones(1,KKX)' * Kx) == (1:1:KKX)' * ones(1,Nx) , 2);
P_y = sum( (ones(1,KKY)' * Ky) == (1:1:KKY)' * ones(1,Ny) , 2);



%Image representation for final synthesis
Liste_pixelsX=create_list(x_,Nx_X,Ny_X);
Liste_pixelsY=create_list(y_,size(y_,1),size(y_,2));


%Mean and Variance for each cluster
[Mean_clusterX Variance_clusterX]=compute_stat_clusters(Liste_pixelsX,reshape(Kx,Nx_X,Ny_X,1),KKX,Nx_X,Ny_X,VARIANCE_COLOR,VARIANCE_SPATIAL);
[Mean_clusterY Variance_clusterY]=compute_stat_clusters(Liste_pixelsY,reshape(Ky,Nx_Y,Ny_Y,1),KKY,Nx_Y,Ny_Y,VARIANCE_COLOR,VARIANCE_SPATIAL);

disp('Creation covariance sparse');

VarianceX=Variance_clusterX{1};
for k=2:KKX,
    VarianceX=blkdiag(VarianceX,Variance_clusterX{k});
end
VarianceX=sparse(VarianceX);

VarianceY=Variance_clusterY{1};
for k=2:KKY,
    VarianceY=blkdiag(VarianceY,Variance_clusterY{k});
end
VarianceY=sparse(VarianceY);

disp('Weight matrix for NN-graph construction');
%Define autocorrelation matrix VX for graph weights between neighbors
% VX_ij is the quadratic Mahalanobis Distance
% d(X_i,X_j)^2 = || X_i - X_j ||^2_E = (X_i - X_j)' E^{-1} (X_i - X_j), where E is the (weighted) empirical covariance matrix
%
% CX2=sum(Mu_x0.^2,1)'*ones(1,KKX);
% VX=CX2+CX2'-2*Mu_x0'*Mu_x0+eye(KKX)*10.;  %we do not want the cluster to be connected with himself when finding near neighbors
%

mean_cluster=(repmat(Mean_clusterX,[1 KKX])-repmat(reshape(Mean_clusterX',1,KKX*5),KKX,1))/ 255.;
weight_cluster=mean_cluster*VarianceX;
weight_cluster=weight_cluster.*mean_cluster;
VX = reshape(sum(reshape(weight_cluster,KKX,5,KKX),2),[KKX KKX]);

Vmax = sort(VX(:),'ascend');
Vmax = Vmax( floor( 90/100 * KKX^2 ) ); % quantile à 90%

% use of a kernel to avoid singularities (also speed-up convergence)
VX( logical(eye(KKX)) ) = inf; %we do not want the cluster to be connected with himself when finding near neighbors
VX = exp(-10*VX/Vmax); % use of quantile 'Vmax' to be robust to the value of parameters 'VARIANCE' and the segmentation



mean_cluster=(repmat(Mean_clusterY,[1 KKY])-repmat(reshape(Mean_clusterY',1,KKY*5),KKY,1))/ 255.;
weight_cluster=mean_cluster*VarianceY;
weight_cluster=weight_cluster.*mean_cluster;
VY = reshape(sum(reshape(weight_cluster,KKY,5,KKY),2),[KKY KKY]);

Vmax = sort(VY(:),'ascend');
Vmax = Vmax( floor( 90/100 * KKY^2 ) ); % quantile à 90%

% use of a kernel to avoid singularities (also speed-up convergence)
VY( logical(eye(KKY)) ) = inf; %we do not want the cluster to be connected with himself when finding near neighbors
VY = exp(-10*VY/Vmax); % use of quantile 'Vmax' to be robust to the value of parameters 'VARIANCE' and the segmentation




%creation of the gradient operator of dimensions nb_voisins x KKX x KKX


GX = create_graph_grad(VX,NB_NEIGHBORS);
GY = create_graph_grad(VY,NB_NEIGHBORS);

grad_X = zeros(KKX); grad_Y = zeros(KKY);
for i=1:NB_NEIGHBORS
    grad_X = grad_X + GX{i};
    grad_Y = grad_Y + GY{i};
end



% create weight on graph to compute average
WX = zeros(KKX); WY = zeros(KKY);
for i=1:NB_NEIGHBORS
    WX = WX + abs(GX{i});
    WY = WY + abs(GY{i});
end
WX = WX ./ repmat(sum(WX,2),[1 KKX]);
WY = WY ./ repmat(sum(WY,2),[1 KKY]);




% Remove spatial coordonates for the clustering used in Optimal Transport
Mu_x = Mu_x(1:3,:);
Mu_y = Mu_y(1:3,:);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2 : Color Matching
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% calcul de la matrice de distance

CX2=sum(Mu_x.^2,1)'*ones(1,KKY);
CY2=ones(KKX,1)*sum(Mu_y.^2,1);
DC_xy=CX2+CY2-2*Mu_x'*Mu_y;
if grey==1
    Mu_y00=Mu_y;
    for ii=1:3
        Mu_y00(ii,:)=0.299*Mu_y(1,:)+0.587*Mu_y(2,:)+0.114*Mu_y(3,:);
    end
    CY2=ones(KKX,1)*sum(Mu_y00.^2,1);
    DC_xy=CX2+CY2-2*Mu_x'*Mu_y00;
end
% dimension des images source/exemple
Nx_X0=size(x0,1);
Ny_X0=size(x0,2);
Nx_Y0=size(y0,1);
Ny_Y0=size(y0,2);

%[Mean_clusterX2 Variance_clusterX2] = compute_stat_clusters(Liste_pixelsX,reshape(Kx,Nx_X,Ny_X,1),KKX,Nx_X,Ny_X,VARIANCE_COLOR,VARIANCE_SPATIAL); % inutile ?
Variance_clusterX2 = Variance_clusterX;
VarianceX2=Variance_clusterX2{1};
for k=2:KKX,
    VarianceX2=blkdiag(VarianceX2,Variance_clusterX2{k});
end

VarianceX2=sparse(VarianceX2);

% precompute data for color transfer synthesis


% average operator from flow for mean and variance color transfer
normalize_rows = @(F) F ./ repmat(max(sum(F ,2),1e-14),[1 size(F,2)]); % so that F * 1 = 1
Mean_flow_x = @(F,Y) normalize_rows(F) *Y;
Mean_flow_y = @(F,X) normalize_rows(F')*X;
Vectorize = @(x) x(:);

% average operator on graph
Mean_graph_x = @(X) WX * X;
Mean_graph_x2 = @(X) WX' * X;
Mean_graph_y = @(Y) WY * Y;
Mean_graph_y2 = @(Y) WY' * Y;

% gradient operator on graph
Grad_graph_x  = @(X) grad_X * X;
Grad_graph_x2 = @(X) grad_X'* X;
Grad_graph_y  = @(Y) grad_Y * Y;
Grad_graph_y2 = @(Y) grad_Y'* Y;



P_x=P_x/Nx;
P_y=P_y/Ny;
normalizeX=diag(P_x(:));
normalizeXm1=normalizeX;
index= normalizeX>0;
normalizeXm1(index)=1./normalizeX(index);

Flow=[];
if USE_SINKHORN
    disp('Sinkhorn Initialization')
    %Initialize with approximated Optimal Transport
    LAMBDA_sinkhorn = 2e2;
    options_sinkhorn.niter = 1e4;
    [cost_sinkhorn,Flow] = sinkhorn_distance(P_x,P_y,DC_xy/max(DC_xy(:)),LAMBDA_sinkhorn,options_sinkhorn);
else
    disp('Optimal Transport initialization')
    %Initialize with true Optimal Transport
    [cost,Flow] = Optimal_Transport(P_x,P_y,DC_xy);
    
end
Flow0=Flow;

DC_xy=DC_xy/max(DC_xy(:));


disp('Optimization')
tic


% projected gradient descent for non-convex problem with the step size paramter 'tau'



% precomputation of normalized matrix Y.Y^T
Mu_yy = Mu_y'*Mu_y / (255.0)^2;


TAU = 3e-3;
%% optimization algorithm
tau = TAU/(1.+2*LAMBDA*max(P_x(:))+2*RHO*(KKX+1)+3*ALPHA*max(P_x(:)));
%tau=1./(1+2*LAMBDA+2*RHO+3*ALPHA);
for it=1:ITMAX,
    
    if do_yreg % matrice de normalisation selon les sommes des colonnes
        normalizeY=diag(P_y);
        
    end
    % CODE for intra-cluster color dispersion (adapted from Nicolas)
    
    
    grad_disp = CY2/255^2 - 2 * normalizeXm1 * Flow * Mu_yy;
    
    
    %     Mx = Mean_graph_x(  normalizeX*( Mean_flow_x(Flow, Mean_clusterY(:,1:3) ) -Mean_clusterX(:,1:3)  )/  255. );
    %
    %
    %     grad_reg=2.*LAMBDA*((Mean_graph_x2(Mx))*Mean_clusterY(:,1:3)')/255.;
    %     if do_yreg
    %         My = Grad_graph_y(normalizeY*( Mean_flow_y(Flow, Mean_clusterX(:,1:3) ) -Mean_clusterY(:,1:3)  )  /  255. );
    %         grad_reg=grad_reg+2.*LAMBDA*(Mean_clusterX(:,1:3)*Grad_graph_y2(My)')/255.;
    %     end
    
    % better regularity term (gradient on graph) : replace mean by difference (gradient)
    
    Mx = Grad_graph_x(normalizeX*( Mean_flow_x(Flow, Mean_clusterY(:,1:3) ) -Mean_clusterX(:,1:3)  )/  255. );
    
    grad_reg=2.*LAMBDA*((Grad_graph_x2(Mx))*Mean_clusterY(:,1:3)')/255.;
    if do_yreg
        My = Grad_graph_y(normalizeY*( Mean_flow_y(Flow, Mean_clusterX(:,1:3) ) -Mean_clusterY(:,1:3)  )  /  255. );
        grad_reg=grad_reg+2.*LAMBDA*(Moyenne_clusterX(:,1:3)*Grad_graph_y2(My)')/255.;
    end
    
    
    grad_fid = 2*RHO * ones(KKX,1)*( sum(Flow,1)./ P_y' );
    
    
    d=DC_xy+grad_disp*ALPHA+grad_reg+grad_fid;
    
    
    
    
    Flow = Flow-d*tau;
    
    for ii=1:KKX,
        Flow(ii,:)=perform_capacity_projection2( Flow(ii,:),ones(1,KKY),P_x(ii));
        
        
    end
    
    
end
toc
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3 : Color Transfer and Image synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% synthetiser les image
disp('Synthesis of images from clusters matchings');
Liste_pixelsX0=create_list(x0,Nx_X0,Ny_X0);
tic
if USE_FAST_SYNTHESIS
    
    XY_ = synthesis_fast(Liste_pixelsX0,Mean_clusterX,VarianceX2,Mean_clusterY,Flow,Nx_X0,Ny_X0,KKX);
else
    XY_ = synthesis(Liste_pixelsX0,Mean_clusterX,Variance_clusterX,Mean_clusterY,Flow,Nx_X0,Ny_X0,KKX);
end
toc
%%
xy_ = reshape(XY_',Nx_X0,Ny_X0,3);



%% fast iterated filter on transport map using guided filter
disp('Post-processing');
tic

addpath('guided_filtering_code/')

N_gf = 10;

I = double(x0) / 255; % guide = original image
p = double(xy_) / 255 - I; % image to be modified = OT color map

r = 1;
eps = 0.1^2;

q = zeros(size(x0));
for k=1:N_gf
 
    q(:, :, 1) = guidedfilter(I(:, :, 1), p(:, :, 1), r, eps);
    q(:, :, 2) = guidedfilter(I(:, :, 2), p(:, :, 2), r, eps);
    q(:, :, 3) = guidedfilter(I(:, :, 3), p(:, :, 3), r, eps);
    
    p = q;
end

xy_ = 255 * (p + I);




toc




%processed image display
figure;
imshow(uint8(xy_)); title(['source->target processing']);drawnow;%








