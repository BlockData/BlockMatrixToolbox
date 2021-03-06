%TESTJACOBIBLOCKPOWER  One-line description here, please.
%
%   output = testJacobiPowerBlock(input)
%
%   Example
%   testJacobiPowerBlock
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-03-02,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.

% initialise le generateur de nombres aleatoires
rng(100);

% create block dimensions
mdims = BlockDimensions({4, [2 3 2]});

% create random data, rounded to avoid floating point errors
data0 = rand(4, 7);
data0 = round(data0 * 1000) / 1000;

% create block matrix instance
data = BlockMatrix.create(data0, mdims);

% display the BlockMatrix
disp('Block matrix:');
disp(data);


%% Cree une bloc-matrice pour les vecteurs
% On utilise une bloc-matrice de 3 blocs, chaque bloc contenant un vecteur.
% La longueur des vecteurs est de (2,3,2).

% create block dimensions
vdims = BlockDimensions({[2 3 2], 1});

% create block matrix instance
tt = BlockMatrix.create(rand(7, 1), vdims);

% display the block-vector (transposed)
disp('Transpose of input vector t:');
disp(tt');


%% Compute problem data

% compute the block-matrix corresponding to maxbet algorithm
AA = blockProduct_uu(data', data);

% ensure symmetry of the matrix...
AA.data = max(AA.data, AA.data');

% create new BlockMatrix representing the normalized input vectors
qq = blockProduct_hs(1./blockNorm(tt), tt);


%% Initialization of AlgoState instance

% creates the algorithm class
algo = JacobiBlockPower(AA);

% create an initial state for the block power iteration algorithm
state0 = BlockPowerAlgoState(qq);

% display content of Algo instance:
disp('Algorithm:');
disp(algo);

% display content of AlgoState instance:
disp('Initial state:');
disp(state0);
disp(state0.vector');


% init residual
resid = 1;
nIters = 10;
residArray = zeros(10, 1);
resArray = cell(nIters, 1);
structArray = [];

tol = 1e-18;


%% Iteration example

stateList = cell(1, nIters);

% iterate until residual is acceptable
state = state0;
for iIter = 1:nIters
    % performs one iteration, and get residual
    state = algo.iterate(state);
    stateList{iIter} = state;
end

% display content of AlgoState instance:
disp('Final state after N iterations:');
disp(state);
disp(state.vector');


%% Solve problem using "solve" method

% iterate until a stopping criterium is met
state = solve(algo, state0);

% display content of AlgoState instance:
disp('Final state using solve:');
disp(state);
disp(state.vector');


%% Compute the path to the solution

% compute the list of positions leading to solution
pathToSolution = convergence(algo, state0);

% display graph of residuals
residuals = cellfun(@(x) x.residual, pathToSolution);
figure; set(gca, 'fontsize', 14);
plot(residuals);
title('Residuals');
xlabel('iterations');

% display the final solution
solution = pathToSolution{end}.vector;
disp('Transposed solution:');
disp(solution');
