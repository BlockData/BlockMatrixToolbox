function X = blockProduct_ss(lambda,A)
% compute 'ss'-type product for block matrices.
%
% AA = blockProduct_ss(lambda, A)
% It corresponds to scalar product along blocks, and scalar product within
% blocks. 
%
% lambda must be a one Block matrix containing a scalar Block.
%
% Example
%   A = BlockMatrix.create(reshape(1:16, [8 2]), {[4 4], [1 1]});
%   disp(A);
%   lambda = BlockMatrix.oneBlock(3);
%   AA = blockProduct_ss(lambda, A)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
[nL, pL] = size(getMatrix(lambda));
kL = blockSize(lambda,1);
kbarL = blockSize(lambda,2);

if or(all([nL,pL]~=[1 1]), all([kL,kbarL]~=[1 1]))
    error(' lambda must be a scalar as one block-Matrix');
end

% create X
dimsX = blockDimensions(A);
matA = getMatrix(A);
scalarlambda = getMatrix(lambda);
X = BlockMatrix.create(scalarlambda*matA, dimsX);