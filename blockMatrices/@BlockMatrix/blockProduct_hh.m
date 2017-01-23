function X = blockProduct_hh(A, B)
% Compute 'hh'-type product for block matrices.
%
% It corresponds to Hadamard product along blocks, and Hadamard product
% within blocks. 
%
% Block-matrices A and B must have the same block-dimension.
%
% Example
%     A = BlockMatrix.create(reshape(1:36, [6 6]), {[3 3], [2 2 2]});
%     B = BlockMatrix.create(magic(6), {[3 3], [2 2 2]});
%     X = blockProduct_hh(A, B)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
if all(blockDimensions(A) ~= blockDimensions(B))
    error('Block dimensions of block-matrices A and B must be the same');
end

% create X
matA = getMatrix(A);
matB = getMatrix(B);
matX = matA .* matB;
X = BlockMatrix.create(matX, blockDimensions(A));
