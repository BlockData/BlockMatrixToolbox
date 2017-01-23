function X = blockProduct_hk(A, B)
% compute 'hk'-type product for block matrices.
%
% It corresponds to Hadamard product along the blocks and Kronecker product
% within the blocks. 
%
% A and B must have the same block-Size.
%
% Example:
%     A = BlockMatrix.create(reshape(1:12, [3 4]), 3, [1 1 1 1]);
%     disp(A);
%     X = blockProduct_hk(A,A)
%
% Reference: A. R. Horn, R. Mathias. Linear Algebra and its Applications,
% 172, pp.337-346 (1992) 
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions
if any(blockSize(A) ~= blockSize(B))
    error('BlockSize of A and BlockSize of B must be the same');
end

% allocate memory for X
prowsX = blockDimensions(A,1) .* blockDimensions(B,1);
pcolX = blockDimensions(A,2) .* blockDimensions(B,2);
X = BlockMatrix.zeros(BlockDimensions({prowsX, pcolX}));

% initialize values for X
for iBlock = 1:blockSize(A, 1)
    for jBlock = 1:blockSize(A, 2)
        BlockA = getBlock(A, iBlock, jBlock);
        BlockB = getBlock(B, iBlock, jBlock);
        BlockX = kron(BlockA, BlockB);
        
        % assign result
        setBlock(X, iBlock, jBlock, BlockX);
    end
end

