function X = blockProduct_sh(A, B)
% compute 'sh'-type product for block matrices.
%
% It corresponds to scalar product along the blocks and Hadamard
% product within the blocks.
% It correponds to multplying a oneBlock Matrix A by a block matrix B
%
% Block-dimension of B must be uniform with the size of each block equals
% size(A).
%
% Example:
%   A = BlockMatrix.oneBlock(2*ones(2,3));
%   B = BlockMatrix.create(reshape(1:36, [4 9]), {[2 2], [3 3 3]});
%   X = blockProduct_sh(A,B)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
kB = blockSize(B,1);
kbarB = blockSize(B,2);
p1 = IntegerPartition(size(A, 1) * ones(1,kB));
p2 = IntegerPartition(size(A, 2) * ones(1,kbarB));
Bdimtest = BlockDimensions({p1, p2});

if all(Bdimtest ~= blockDimensions(B))
    error('Block-dimension of B must be uniform with the size of each block equals size(A)');
end

% create X
matA = getMatrix(A);
X = BlockMatrix.zeros(blockDimensions(B));

% Compute Blocks of X
for iBlock = 1:1:blockSize(B, 1)
    for jBlock = 1:blockSize(B, 2)
        % extract blocks of the two input block matrices
        blockB = getBlock(B, iBlock, jBlock);
        % compute 'h'-product of A and blocks of B
        XBlock = matA .* blockB;
        % assign result
        setBlock(X, iBlock, jBlock, XBlock);
    end
end
