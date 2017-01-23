function X = blockProduct_su(A, B)
% compute 'su'-type product for block matrices.
%
% It corresponds to scalar product along the blocks, and usual product
% within the blocks. 
% It correponds to multiplying a oneBlock-Matrix A by a block matrix B
% rowblocks of B must be uniform with each parts equals columns of A.
%
% Example:
%   A = BlockMatrix.oneBlock(ones(2,3));
%   B = BlockMatrix.create(reshape(1:36, [6 6]), {[3 3], [2 2 2]});
%   disp(B);
%   X = blockProduct_su(A, B);
%
% Reference: M. Günther, L. Klotz. Linear Algebra and its Applications.
% 437, pp. 948-956,(2012) 
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
kB = blockSize(B,1);
browstest = IntegerPartition(size(getMatrix(A),2)*ones(1,kB));
browsB = IntegerPartition(blockDimensions(B,1));
if or(all(blockSize(A)~=[1, 1]), all(browstest~=browsB))
    error('row blocks of B must be uniform with each part equals size(A,2)');
end

% create X
matA = getMatrix(A);
browsX = IntegerPartition(size(A,1) * ones(1, kB));
bcolX = blockDimensions(B,2);
X = BlockMatrix.zeros(BlockDimensions({browsX, bcolX}));

% Compute Blocks of X
for iBlock = 1:blockSize(B, 1)
    for jBlock = 1:1:blockSize(B, 2)
        % extract blocks of block matrix B
        blockB = getBlock(B, iBlock, jBlock);
        % compute 'u'-product of A and blocks of blockB
        XBlock = matA * blockB;
        % assign result
        setBlock(X, iBlock, jBlock, XBlock);
    end
end
