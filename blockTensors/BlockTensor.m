classdef BlockTensor < handle
%BLOCKTENSOR Parent interface for block-partitioned multidimensional tensors
%
%   See also
%     BlockMatrix, BlockDiagonal
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-11-03,    using Matlab 9.1.0.441655 (R2016b)
% Copyright 2016 INRA - BIA-BIBS.


% no properties are defined, nor constructor

%% Interface Methods
% Methods in this cell are declared for implementation in sub-classes

methods (Abstract)
    % Return the dimensions of the block in the specified dimension
    dims = blockDimensions(this, dim)
    
    % Return the number of dimensions of this block tensor
    dim = dimensionality(this)
       
    % Return the total number of blocks in this block tensor
    n = blockNumber(this, varargin)
    
    % Return the number of blocks in each dimension
    n = blockSize(this, varargin)
    
end % end abstract methods


end % end classdef

