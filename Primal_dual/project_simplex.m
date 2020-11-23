 function p = project_simplex(x, eta, dir)
%function p = project_simplex(x, eta, dir)
%
% This procedure computes the projection onto the constraint set:
%
%                  x => 0   AND   1'x = eta
%
% When the input 'x' is an array, the computation can vary as follows:
%  - dir = 0 --> 'x' is processed as a single vector [DEFAULT]
%  - dir > 0 --> 'x' is processed block-wise along the specified direction
%
%  INPUTS
% ========
%  x   - ND array
%  eta - positive, scalar or ND array compatible with the blocks of 'x'
%  dir - integer, direction of block-wise processing
% 
%  DEPENDENCIES
% ==============
%  prox_max.m  - located in the folder 'multi'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version : 1.0 (27-04-2017)
% Author  : Giovanni Chierchia
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2017
%
% This file is part of the codes provided at http://proximity-operator.net
%
% By downloading and/or using any of these files, you implicitly agree to 
% all the terms of the license CeCill-B (available online).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% default inputs
if nargin < 3 || (~isempty(dir) && dir == 0)
    dir = [];
end

% check input
sz = size(x); sz(dir) = 1;
if any( eta(:) <= 0 ) || ~isscalar(eta) && any(size(eta) ~= sz)
    error('''eta'' must be positive and either scalar or compatible with the blocks of ''x''')
end
%-----%


% compute the projection
p = x - prox_max(x, eta, dir);