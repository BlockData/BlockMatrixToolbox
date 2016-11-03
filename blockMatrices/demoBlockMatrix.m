%DEMOBLOCKMATRIX Simple demo file for manipulation of BlockMatrix objects
%
%   Example
%   demoBlockMatrix
%
%   See also
%     BlockMatrix
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-02-19,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - Cepia Software Platform.

% create data 
data = reshape(1:28, [4 7]);

% specifies dimensions of the block matrix:
% 2 rows (with 2 and 2 rows)
% 3 columns (with 2, 3, and 2 columns resepctively)
parts = {[2 2], [2 3 2]};

% create the BlockMatrix object
BM = BlockMatrix(data, parts);

%% display info about the block Matrix
disp(BM);

%% Create a new block matrix with different size
dims2 = {[1 2 1], [4 3]};
BM2 = BlockMatrix(BM, dims2);
disp(BM2);

%% show the content of one of the blocks
getBlock(BM, 2, 2)

%% Also display the content of all blocks
displayBlocks(BM)