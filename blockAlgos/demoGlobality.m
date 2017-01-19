function varargout = demoGlobality(varargin)
%DEMOGLOBALITY  One-line description here, please.
%
%   output = demoGlobality(input)
%
%   Example
%   demoGlobality
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2017-01-19,    using Matlab 9.1.0.441655 (R2016b)
% Copyright 2017 INRA - Cepia Software Platform.


%% Initialisations

% initialise le generateur de nombres aleatoires
rng(100);

% create block dimensions
mdims = BlockDimensions({4, [2 3 2]});


%% Initialisations des strucutres de donnees

% create random data, rounded to avoid floating point errors
data0 = rand(4, 7);
data0 = round(data0 * 1000) / 1000;

% create block matrix instance
data = BlockMatrix.create(data0, mdims);

% display the BlockMatrix
disp('Block matrix:');
disp(data);

% compute the block-matrix corresponding to maxbet algorithm
AA = blockProduct_uu(data', data);

% ensure symmetry of the matrix...
AA.data = max(AA.data, AA.data');

% creates the algorithm class
algo = JacobiBlockPower(AA);


%% Globalite simple

res = globality(algo)


%% Globalite simple

res = globality(algo)

