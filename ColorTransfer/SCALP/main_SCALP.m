% This code is free to use for any non-commercial purposes.
% It contains an implementation of the superpixel method proposed in:
% [1] - Rémi Giraud, Vinh-Thong Ta and Nicolas Papadakis
%       SCALP: Superpixels with Contour Adherence using Linear Path
%       International Conference on Pattern Recognition (ICPR), 2016
%
% If you use this code, please cite [1].
%
% A contour prior map can be added to our method. 
% The contour detection used method in [1] is available with guidelines at 
% https://github.com/pdollar/edges
% Other contour detection methods can be found here
% https://www2.eecs.berkeley.edu/Research/Projects/CS/vision/grouping/resources.html
%
% (C) Rémi Giraud, 2017
% rgiraud@labri.fr, www.labri.fr/~rgiraud/downloads
% University of Bordeaux, LaBRI, IMB


% clear all
% close all
% clc

%Compilation
mex -O CFLAGS="\$CFLAGS -Wall -Wextra -W -std=c99" SCALP_v1.c -outdir ./
%%

% Image loading
I = imread('../images/parrot-1.jpg');
[h,w,z] = size(I);

% Parameters
K = 300;    % Superpixel number
m = 0.01;  % Compactness parameter

% Contour detection
% Plug here your contour detector C (see readme.txt file)
% contour detector must be between 0 and 1
% C = contour_detection...

% SCALP
S = SCALP_v1(single(I)/255,K,m);
% S = SCALP_v1(single(I)/255,K,m,single(C));  %with contour prior

% Plotting results
B = zeros(h,w);
for i=1:max(S(:))
    j = (S == i);
    bb = bwboundaries(j);
    if ~isempty(bb)
        for k=1:length(bb{1})
            B(bb{1}(k,1),bb{1}(k,2)) = 1;
        end
    end
end

figure,
imagesc(double(I)/255.*repmat(~B,[1 1 3])); 



