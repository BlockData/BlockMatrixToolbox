function X = blockProduct_uk(A, B)
% compute 'uk'-typeproduct for block matrices.
%
% It corresponds to usual product along the blocks and Kronecker
% product within the blocks.
%
% Block-rows of B and block-columns of A must be uniform.
%
% Example:
%   A = BlockMatrix.create(reshape(1:12, [3 4]), {[3], [1 1 1 1 ]});
%   B = BlockMatrix.create(reshape(1:12, [4 3]), {[1 1 1 1], [1 2]});
%   X = blockProduct_uk(A,A')
%
% References: 	W. De Launey J. Seberrya Journal of Combinatorial Theory,
% Series A. 66 (2), pp. 192–213 (1994).

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on block dimensions
pA = size(getMatrix(A),2);
nB = size(getMatrix(B),1);
kB = blockSize(B,1);
kbarA = blockSize(A,2);
kbarB = blockSize(B,2);

s1 = or(round((nB/kB))~=(nB/kB),round((pA/kB)~=(pA/kB)));
p1 = IntegerPartition(round((nB/kB))*ones(1,kB));
p2 = IntegerPartition(round((pA/kB))*ones(1,kbarA));
s2 = or(p1~=blockDimensions(B,1),p1~=blockDimensions(A,2));
if not(and(isUniform(blockDimensions(A,2)),isUniform(blockDimensions(B,1))))
    error(' block rows of B and block columns of A must be uniform');
end

% create X
prowsX = term(blockDimensions(B,1),1)*blockDimensions(A,1);
pcolX = term(blockDimensions(A,2),1)*blockDimensions(B,2);
X = BlockMatrix.zeros(BlockDimensions({prowsX, pcolX}));

for iBlock = 1:blockSize(A, 1)
    for jBlock = 1:blockSize(B, 2)
        % compute (i,j) blocks
        s = zeros(prowsX(iBlock),pcolX(jBlock));
        for kBlock = 1:blockSize(A, 2)
            % extract blocks of the two input block matrices
            BlockA = getBlock(A, iBlock, kBlock);
            BlockB = getBlock(B, kBlock, jBlock);
            s = s + kron(BlockA, BlockB);
        end
        % assign result
        setBlock(X, iBlock, jBlock, s);
    end
end
