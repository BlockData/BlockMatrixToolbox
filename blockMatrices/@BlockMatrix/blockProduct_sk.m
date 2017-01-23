function X = blockProduct_sk(A, B)
% compute 'sk'-type product for block matrices.
%
% It corresponds to scalar product along blocks, and kronecker product
% within blocks. 
%
% A must be oneBlock-Matrix.
%
% Example
%     B = BlockMatrix.create(reshape(1:36, [6 6]), {[3 3], [2 2 2]})
%     A = oneBlock(reshape(1:4,[2,2]))
%     X = blockProduct_sk(A,B);
%
% Reference: R. H. Koning, H. Neudecker†, T.  Wansbeek. Linear Algebra and
% its Appl. 149, pp 165–184,(1991) 

% ------
% Author: Mohamed Hanafi
% e-mail: mohamed.hanafi@oniris-nantes.fr
% Created: 2015-06-00,  using Matlab(R2015a)

% check conditions on dimensions

if  all(blockSize(A)~=[1,1])
    error('A must be one Block-Matrix');
end

browsX = IntegerPartition(size(getMatrix(A),1)*blockDimensions(B,1));
bcolX = IntegerPartition(size(getMatrix(A),2)*blockDimensions(B,2));

% create X
matA = getMatrix(A);
X = BlockMatrix.zeros(BlockDimensions({browsX, bcolX}));

% Compute Blocks of X
for iBlock = 1:blockSize(B, 1)
    for jBlock = 1:blockSize(B, 2)
        % extract blocks of block matrix B
        blockB = getBlock(B, iBlock, jBlock);
        % compute 'k'-product of A and blocks of B
        XBlock = kron(matA, blockB);
        % assign result
        setBlock(X, iBlock, jBlock, XBlock);
    end
end
