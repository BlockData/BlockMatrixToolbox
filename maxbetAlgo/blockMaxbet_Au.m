function u = blockMaxbet_Au(A, u0)
%BLOCKMAXBET_AU Maxbet algorithm: computes vector that maximizes u'*A*u
%
%   U = maxbet_Au(A, U0)
%   A is a square blockmatrix representing X' * X
%   U0 is a block vector representing the initialisation vector.
%   The function computes the vector that maximizes u'*A*u subject to
%   ||u_k|| = 1. 
%
%   Example
%   
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-06-13,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.


%% ensure block norm of vector

u = u0;

% create new BlockMatrix representing the normalized input vector
u = blockProduct_hs(1 ./ blockNorm(u), u);


%% iterate
for iIter = 1:100
    % compute block-product
    v = blockProduct_uu(A, u);
    
    % block normalization of vector
    v = blockProduct_hs(1 ./ blockNorm(v), v); 
    
%     % compute residuals
%     resid = norm(blockNorm(u) - blockNorm(v));
    
    u = v;
end
