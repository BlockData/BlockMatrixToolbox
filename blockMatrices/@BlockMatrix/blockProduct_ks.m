function X = blockProduct_ks(A, B)
% compute 'ks'-type product for block matrices.
%
% It corresponds to kronecker product along the blocks and scalar product
% within the blocks. 
% 
% A must be a scalar Block Matrix and B must be one Block-Matrix.
%
% Example
%   A = BlockMatrix.scalarBlock(reshape(1:4, [2 2]));
%   disp(A);
%   B = BlockMatrix.oneBlock(reshape(1:4, [2 2]));
%   X = blockProduct_ks(A,B)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
[nA, pA] = size(A);
kA      = blockSize(A,1);
kbarA   = blockSize(A,2);
kB      = blockSize(B,1);
kbarB   = blockSize(B,2);

if any([nA,pA] ~= [kA,kbarA])
    error('Blocks of A must be a scalars');
end
if any([kB,kbarB] ~= [1,1])
    error('B must be one block-matrix');
end

% create X
[nB,pB] = size(getMatrix(B));
prows = IntegerPartition(nB*ones(nA,1));
pcol = IntegerPartition(pB*ones(pA,1));

matA = getMatrix(A);
matB = getMatrix(B);
X = BlockMatrix.create(kron(matA,matB), BlockDimensions({prows, pcol}));
