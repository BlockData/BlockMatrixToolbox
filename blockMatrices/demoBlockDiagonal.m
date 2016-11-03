%DEMOBLOCKDIAGONAL  Simple demo file for manipulation of BlockDiagonal objects
%
%   output = demoBlockDiagonal(input)
%
%   Example
%   demoBlockDiagonal
%
%   See also
%     BlockDiagonal
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-03-10,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - Cepia Software Platform.


%% Creation d'un objet BlockDiagonal

% create two rectangular matrices
mat1 = reshape(1:6, [2 3]);
mat2 = reshape(7:12, [3 2]);

% use the blkdiag to show an example of block diagonal stored as a 2D array
blkdiag(mat1, mat2)

BD = BlockDiagonal(mat1, mat2);
disp(BD);


%% Affiche les blocs

displayBlocks(BD)