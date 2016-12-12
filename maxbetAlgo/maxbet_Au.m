function u = maxbet_Au(A, u0, vdims)
%MAXBET_AU Maxbet algorithm: computes vector that maximizes u'*A*u subject to |u_k|=1
%
%   U = maxbet_Au(A, U0, VDIMS)
%
%   Example
%   maxbet_Au
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-06-13,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.


%% pre-compute some data

nBlocks = length(vdims);
blockInds = cell(1, nBlocks);

inds = 1:length(u0);

pos1 = 1;
for i = 1:nBlocks
    pos2 =  pos1 + vdims(i) - 1;
    blockInds{i} = inds(pos1:pos2);
    pos1 = pos2 + 1;
end


%% ensure block norm of vector

u = u0;

for iBlock = 1:nBlocks
    inds = blockInds{iBlock};
    uk = u(inds);
    ukn = uk / norm(uk, 2);
    u(inds) = ukn;
end


%% iterate

v = zeros(size(u));

for iIter = 1:100
    % compute v_k = A * u_k
    for iBlock = 1:nBlocks
        inds = blockInds{iBlock};
        
        % compute v_i = sum_j A_ij * u_j
        vk = zeros(size(u(inds)));
        for jBlock = 1:nBlocks
            inds2 = blockInds{jBlock};
            vk = vk + A(inds, inds2) * u(inds2);
        end
        v(inds) = vk;
    end
    
    % rescale each v_k
    for iBlock = 1:nBlocks
        inds = blockInds{iBlock};
        vk = v(inds);
        vkn = vk / norm(vk, 2);
        v(inds) = vkn;
    end
    
    % assign u = v
    u = v;
end
