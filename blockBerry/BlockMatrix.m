classdef BlockMatrix < handle
%BLOCKMATRIX Matrix that can be divided into several blocks
%
%   Class BlockMatrix
%
%   Example
%   BlockMatrix
%
%   See also
%     BlockDimension

% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2015-02-19,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % contains the dimensions of blocs in each dimension, as a cell array
    % containing row vectors of integers
    dims;
    
    % contains the data, as a vector or an array. 
    % Idea is to access only via linear indexing
    data;
    
end % end properties


%% Constructor
methods
    function this = BlockMatrix(varargin)
        % Constructor for BlockMatrix class
        %
        %   data = reshape(1:28, [4 7]);
        %   BM = BlockMatrix(data, {[2 2], [2 3 2]});
        %
        
        if isempty(varargin)
            % populate with default data
            this.data = 1:28;
            this.dims = {[2 2], [2 3 2]};
            
        elseif nargin == 2
            if isnumeric(varargin{1})
                % initialsation constructor
                this.data = varargin{1};
                
            elseif isa(varargin{1}, 'BlockMatrix')
                % copy constructor
                bm = varargin{1};
                this.data = bm.data;
            else
                error('first argument must be a matrix or a block matrix');
            end
            
            % second argument represents block dimension. This can be
            % either a BlockDimension object, or a cell array containing
            % the block sizes in each dimension.
            var2 = varargin{2};
            if isa(var2, 'BlockDimensions')
                this.dims = var2;
            elseif iscell(var2)
                this.dims = BlockDimensions(var2);
            else
                error('second argument must be a cell array or a BlockDimensions object');
            end
            
        else
            error('Requires two input arguments');
        end

    end

end % end constructors


%% Methods specific to BlockMatrix object
% sorted approximately from high-level to low-level

methods
    
    function block = getBlock(this, row, col)
        % return the (i-th, j-th) block 
        %
        %   BLK = getBlock(BM, ROW, COL)
        %
        
        % determine row indices of block rows
        parts1 = getBlockDimensions(this.dims, 1);
        rowInds = (1:parts1(row))' + sum(parts1(1:row-1));

        % determine column indices of block columns
        parts2 = getBlockDimensions(this.dims, 2);
        colInds = (1:parts2(col)) + sum(parts2(1:col-1));
        
        % compute full size of block matrix
        dim = [sum(parts1) sum(parts2)];
        
        % compute indices of block elements in data
        colInds2 = repmat(colInds, length(rowInds), 1);
        rowInds2 = repmat(rowInds, 1, length(colInds));
        inds = sub2ind(dim, rowInds2, colInds2);
        
        % extract data element corresponding to block. 
        block = this.data(inds);
    end
    
    function setBlock(this, row, col, blockData)
        % set the data for the (i-th, j-th) block 
        %
        %   setBlock(BM, ROW, COL, DATA)
        %
        
        % determine row indices of block rows
        parts1 = getBlockDimensions(this.dims, 1);
        rowInds = (1:parts1(row))' + sum(parts1(1:row-1));
        
        % check number of rows of input data
        if length(rowInds) ~= size(blockData, 1)
            error('block data should have %d rows, not %d', ...
                length(rowInds), size(blockData, 1));
        end

        % determine column indices of block columns
        parts2 = getBlockDimensions(this.dims, 2);
        colInds = (1:parts2(col)) + sum(parts2(1:col-1));
        
        % check number of columns of input data
        if length(colInds) ~= size(blockData, 2)
            error('block data should have %d columns, not %d', ...
                length(colInds), size(blockData, 2));
        end

        % compute full size of block matrix
        dim = [sum(parts1) sum(parts2)];
        
        % compute indices of block elements in data
        colInds2 = repmat(colInds, length(rowInds), 1);
        rowInds2 = repmat(rowInds, 1, length(colInds));
        inds = sub2ind(dim, rowInds2, colInds2);
        
        % extract data element corresponding to block. 
        this.data(inds) = blockData;
    end
end

%% Methods that depends uniquely on BlockDimensions object

methods
    function dims = getBlockDimensions(this, dim)
        % Return the dimensions of the block in the specified dimension
        %
        % DIMS = getBlockDimensions(BM, IND)
        %
        dims = getBlockDimensions(this.dims, dim);
    end
    
    function dim = getDimensionality(this)
        % Return the number of dimensions of this block matrix (usually 2)
        dim = getDimensionality(this.dims);
    end
    
    function siz = getSize(this, varargin)
        % Return the size in each direction of this block matrix object
        siz = getSize(this.dims, varargin{:});
    end
    
    function n = getBlockNumber(this, varargin)
        % Return the total number of blocks in this block matrix, or the
        % number of blocks in a given dimension
        %
        % N = getBlockNumber(BM);
        % N = getBlockNumber(BM, DIM);
        %
        n = getBlockNumber(this.dims, varargin{:});
    end
    
    function n = getBlockNumbers(this)
        % Return the number of blocks in each dimension
        n = getBlockNumbers(this.dims);
    end
end


%% Overload some native methods

methods
    function res = transpose(this)
        % transpose this BlockMatrix
        res = ctranspose(this);
    end
    
    function res = ctranspose(this)
        % overload the transpose operator for BlockMatrix object
        
        % ensure the new matrix is a rectangular array, with the new size
        siz = getSize(this);
        data2 = reshape(this.data, siz)';
        
        % transpose the BlockDimensions object, and create nesw BlockMatrix
        dims2 = transpose(this.dims);
        res = BlockMatrix(data2, dims2);
    end
    
    function res = cat(dim, varargin)
        % overload concatenation method for arbitrary dimension (between 1 and 2...) 
        switch dim
            case 1
                res = vertcat(varargin{:});
            case 2
                res = horzcat(varargin{:});
            otherwise
                error('unsupported dimension: %d', dim);
        end
    end
    
    function res = horzcat(this, varargin)
        % Overload the horizontal concatenation operator
        
        % initialize block dimension and data to that of first BlockMatrix
        data2 = reshape(this.data, getSize(this));
        dims2 = this.dims;
        
        for i = 1:length(varargin)
            var = varargin{i};
            
            dataToAdd = reshape(var.data, getSize(var));
            if size(dataToAdd, 1) ~= size(data2, 1)
                error('BlockMatrices should have same number of rows');
            end
            
            data2 = [data2 dataToAdd]; %#ok<AGROW>
            dims2 = [dims2 var.dims]; %#ok<AGROW>
        end
        
        res = BlockMatrix(data2, dims2);
    end
    
    function res = vertcat(this, varargin)
        % Overload the vertical concatenation operator
        
        % initialize block dimension and data to that of first BlockMatrix
        data2 = reshape(this.data, getSize(this));
        dims2 = this.dims;
        
        for i = 1:length(varargin)
            var = varargin{i};
            
            dataToAdd = reshape(var.data, getSize(var));
            if size(dataToAdd, 2) ~= size(data2, 2)
                error('BlockMatrices should have same number of columns');
            end
            
            data2 = [data2 ; dataToAdd]; %#ok<AGROW>
            dims2 = [dims2 ; var.dims]; %#ok<AGROW>
        end
        
        res = BlockMatrix(data2, dims2);
    end
    
end

%% Display methods

methods
    
    function disp(this)
        % display the content of this BlockMatrix object
        
        % loose format: display more empty lines
        isLoose = strcmp(get(0, 'FormatSpacing'), 'loose');
        
        % get BlockMatrix total size
        dim = getSize(this);
        nRows = dim(1);
        nCols = dim(2);
        
        % Display information on block sizes
        disp(sprintf('BlockMatrix object with %d rows and %d columns', nRows, nCols)); %#ok<DSPS>
        parts1 = getBlockDimensions(this.dims, 1);
        disp(sprintf('  row dims: %s', formatParts(parts1))); %#ok<DSPS>
        parts2 = getBlockDimensions(this.dims, 2);
        disp(sprintf('  col dims: %s', formatParts(parts2))); %#ok<DSPS>
        
        if isLoose
            fprintf('\n');
        end
        
        function string = formatParts(parts)
            % display parts as a list of ints, or as empty
            if isempty(parts)
                string = '(empty)';
            else
                pattern = strtrim(repmat(' %d', 1, length(parts)));
                string = sprintf(pattern, parts);
            end
        end
    end
    
    
    function displayBlocks(this)
        % get BlockMatrix total size
        nRowBlocks = length(getBlockDimensions(this.dims, 1));
        nColBlocks = length(getBlockDimensions(this.dims, 2));
        
        for row = 1:nRowBlocks
            for col = 1:nColBlocks
                disp(sprintf('Block (%d,%d)', row, col)); %#ok<DSPS>
                disp(getBlock(this, row, col));
            end
        end
        
    end
    
end % end methods

end % end classdef

