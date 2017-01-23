function X = blockProduct_ku(A, B)
% compute 'ku'-type product for block matrices.
%
% It corresponds to kronecker product along the blocks, and usual product
% within the blocks. 
%
% Transposition rule: blockProduct_ku(A, B)=(blockProduct_su(B', A'))'
%
% Example:
%   A = BlockMatrix.create(reshape(1:36, [6 6]), {[2 2 2], [3 3]});
%   B = BlockMatrix.oneBlock(ones(3,2));
%   X = blockProduct_ku(A,B);
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
AA = B';
BB = A';
kBB = blockSize(BB,1);
browstest = IntegerPartition.ones(kBB) * size(AA,2);
browsBB = IntegerPartition(blockDimensions(BB,1));
if any(blockSize(AA) ~= [1 1]) || any(browstest ~= browsBB)
    error('column blocks of A must be uniform with each part equals size(B,1)');
end

X = blockProduct_su(AA, BB);
X = X';
