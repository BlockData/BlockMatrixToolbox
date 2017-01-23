function X = blockProduct_hs(A, B)
% compute 'hs'-type product for block matrices.
%
% It corresponds to hadamard product along the blocks and scalar product
% within the blocks. 
%
% A must be a scalarBlock-Matrix and the block Size of A and B must be
% same. 
%
% Example:
%    A = BlockMatrix.scalarBlock(reshape(1:4, [2 2]))
%    B = BlockMatrix.create(reshape(1:16, [4 4]), {[2 2], [2 2]});
%    X = blockProduct_hs(A,B)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on block dimensions
[nA, pA] = size(getMatrix(A));
kA = blockSize(A,1);
kbarA = blockSize(A,2);
kB = blockSize(B,1);
kbarB = blockSize(B,2);

if or(all([kA kbarA]~=[kB kbarB]),all([kA kbarA]~=[nA pA]))
    error('A must be a scalar Block matrix and the block structure of A and B must be the same');
end

% create blockdims of X
brows = blockDimensions(B,1);
bcol = blockDimensions(B,2);
X = BlockMatrix.zeros(BlockDimensions({brows, bcol}));

% Compute blocks of X
for iBlock = 1:blockSize(B, 1)
    for jBlock = 1:blockSize(B, 2)
        % extract blocks of the two input block matrices
        BlockA = getBlock(A, iBlock, jBlock);
        BlockB = getBlock(B, iBlock, jBlock);
        % compute 'u'-product of blocks
        XBlock = BlockA * BlockB;
        % assign result
        setBlock(X, iBlock, jBlock,XBlock);
    end
end
