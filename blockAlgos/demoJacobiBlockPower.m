%DEMOJACOBIBLOCKPOWER Demonstration of Jacobi Block-Power Iteration algorithm
%
%   output = demoJacobiBlockPower(input)
%
%   Example
%   demoJacobiBlockPower
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-01-28,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.

%% create a block matrix for initial matrix

% create block dimensions
mdims = BlockDimensions({4, [2 3 2]});

% create block matrix instance
data = BlockMatrix(rand(4, 7), mdims);

% display the BlockMatrix
disp('Block matrix:');
disp(data);


%% create a block matrix for vector

% create block dimensions
vdims = BlockDimensions({[2 3 2], 1});

% create block matrix instance
q0 = BlockMatrix(rand(7, 1), vdims);

% display the block-vector (transposed)
disp('Transpose of input vector t:');
disp(q0');


%% Initialize data for block iteration

% compute the block-matrix corresponding to maxbet algorithm
AA = blockProduct_uu(data', data);
reveal(AA);

% create new BlockMatrix representing the normalized input vectors
qq = blockProduct_hs(1./blockNorm(q0), q0);

% display the block-vector (transposed)
disp('initial vector:');
disp(qq');


%% iterate for a given number of iterations


% create the algorithm iterator
algo = JacobiBlockPower(AA, qq);

% allocate memory for result
nIter = 20;
residList = zeros(nIter, 1);
lambdaList = zeros(nIter, 1);

% iterate
for iIter = 1:nIter
    [q, residList(iIter)] = iterate(algo);
    lambdaList(iIter) = eigenValue(algo);
end


% display result
figure; set(gca, 'fontsize', 14); hold on;
plot([1 nIter], lambdaList([end end]), 'k');
plot(lambdaList, 'color', 'b', 'linewidth', 2);
xlabel('Iteration Number');
ylabel('Eigen value estimation');


% % init residual
% resid = 1;
% 
% % iterate until residual is acceptable
% nIter = 0;
% while true
%     % increment iteration
%     nIter = nIter + 1;
%     
%     % performs one iteration, and get residual
%     [q, resid] = iterator.iterate();
%     
%     % check stopping conditions
%     if resid < this.tolerance
%         break;
%     end
%     if nIter >= this.nbMaxIter
%         break;
%     end
% end
