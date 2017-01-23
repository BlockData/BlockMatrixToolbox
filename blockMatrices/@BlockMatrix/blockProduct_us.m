function X = blockProduct_us(A, B)
% compute 'us'-type product for block matrices.
%
% It corresponds to usual product along the blocks and scalar product
% within the blocks. 
%
% Validity conditions:
% * A must be a scalar block-matrix
% * the number of column-blocks of A must be the same as the number of
%   row-blocks of B  
% * the row partition of B must be uniform
%
% Example:
%   A = BlockMatrix.scalarBlock([1 2 3;3 2 1]);
%   B = BlockMatrix.create(reshape(1:36, [4 9])', {[3 3 3], [2 2]});
%   disp(B);
%   X = blockProduct_us(A, B)
%

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on block dimensions
[nA, pA] = size(A);
nB = size(B, 1);
kA = blockSize(A, 1);
kbarA = blockSize(A, 2);
if any([kA kbarA] ~= size(A))
    error('Each block from A must be a scalar');
end
if round((nB/pA)) ~= (nB/pA)
    error('The block structure of A and B must be same');
end
if any(IntegerPartition(round((nB/pA)*ones(1,pA))) ~= blockDimensions(B,1))
    error('The rows of the blocks from B must be of the same size');
end

% create blockdims of res
bdimX = BlockDimensions({IntegerPartition(round((nB/pA))*ones(1,nA)), blockDimensions(B,2)});
X = BlockMatrix.zeros(bdimX);
brows = blockDimensions(B,1);
bcol = blockDimensions(B,2);

for iBlock = 1:blockSize(A, 1)
    for jBlock = 1:blockSize(B, 2)
        % compute (i,j) blocks
        s = zeros(brows(iBlock), bcol(jBlock));
        for kBlock = 1:blockSize(A, 2)
            % extract blocks of the two input block matrices
            BlockA = getBlock(A, iBlock, kBlock);
            BlockB = getBlock(B, kBlock, jBlock);
            s = s + BlockA * BlockB;
        end
        % assign result
        setBlock(X, iBlock, jBlock, s);
    end    
end
