function X = blockProduct_uu(A, B)
% compute 'uu'-type product for block matrices.
%
% It corresponds to usual product along the blocks and usual product within
% the blocks. 
%
% the block-rows of B and block-columns of A  must be equal check
% conditions on dimensions.
%
% Example
%   A = BlockMatrix(reshape(1:36, [6 6]), {[2 2 2], [3 3]});
%   B = BlockMatrix(reshape(1:36, [6 6]), {[3 3], [2 1 2 1]});
%   X = blockProduct_uu(A,B)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
if blockDimensions(A, 2)~= blockDimensions(B, 1)
    error(' block rows of B and block columns of A must be equal');
end

% create X
dimsX = BlockDimensions({blockDimensions(A, 1), blockDimensions(B, 2)});
matA = getMatrix(A);
matB = getMatrix(B);
X = BlockMatrix.create(matA * matB, dimsX);
