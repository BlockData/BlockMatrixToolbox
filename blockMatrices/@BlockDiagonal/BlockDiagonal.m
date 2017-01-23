classdef (InferiorClasses = {?BlockMatrix}) BlockDiagonal <  BlockMatrix
%BLOCKDIAGONAL Block Matrix with zeros blocks except on diagonal blocks
%
%   BlockDiagonal objects are constructed from the list of blocks located
%   on the diagonal. The Block-Dimensions of the block-diagonal is computed
%   automatically.
%
%   Example
%   % create a block diagonal matrix
%   BD = BlockDiagonal({rand(2, 3), rand(2, 2), rand(1, 2)});
%   % transpose and multiplies two blockdiagonals
%   BD' * BD
%
%   See also
%     BlockMatrix, BlockDimensions
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-03-02,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % the set of diagonal blocks, as a cell array containing matrices
    diags;
    
    % the block dimensions of this diagonal matrix, as a cell array of
    % integer partitions.
    dims;
    
end % end properties

%% Static methods
methods (Static)
    function BM = scalarBlock(vect)
        %SCALARBLOCK  Converts a vector to a BlockDiagonal with only scalar blocks
        %
        %   BDM = scalarBlock(VECT)
        %
        %   Example
        %   BDM = BlockDiagonal.scalarBlock([1 2 3]);
        %   reveal(BDM)
        %          1  1  1
        %       1  +  +  +
        %       1  +  +  +
        %       1  +  +  +
        %
        %   See also
        %     BlockDiagonal, BlockMatrix.scalarBlock
        
        BM = BlockDiagonal(num2cell(vect));
    end
end


%% Constructor
methods
    function this = BlockDiagonal(varargin)
        % Constructor for BlockDiagonal class
        %
        %   diagBlocks = {rand(2,3), rand(2,2), rand(3, 2)};
        %   BD = BlockDiagonal(diagBlocks);
        %
        
        if nargin == 0
            error('Require at least one input argument');
        end
        
        if iscell(varargin{1})
            % blocks are given as a cell array of matrices
            this.diags = varargin{1};
            computeDimensions();
            
        elseif all(cellfun(@isnumeric, varargin))
            % blocks are given as varargin
            this.diags = varargin;
            computeDimensions();
            
        elseif isa(varargin{1}, 'BlockDiagonal')
            % copy constructor
            var1 = varargin{1};
            this.diags = var1.diags;
            this.dims = var1.dims;
            
        elseif isa(varargin{1}, 'BlockMatrix')
            % copy constructor from another type of BlockMatrix object
            bm = varargin{1};
            blocks = cell(1, blockNumber(bm));
            iBlock = 0;
            for i = 1:blockSize(bm, 1)
                for j = 1:blockSize(bm, 2)
                    iBlock = iBlock + 1;
                    blocks{iBlock} = getBlock(bm, i, j);
                end
            end
            
            this.diags = blocks;
            computeDimensions();

        else
            error('input argument must be a cell array of matrices');
        end
        
        
        function computeDimensions()
            % inner function that computes block dimensions
            nDiags = length(this.diags);
            dims1 = zeros(1, nDiags);
            dims2 = zeros(1, nDiags);
            for iDiag = 1:nDiags
                siz = size(this.diags{iDiag});
                dims1(iDiag) = siz(1);
                dims2(iDiag) = siz(2);
            end
            this.dims = BlockDimensions({dims1, dims2});
        end
    end

end % end constructors


%% Methods specific to BlockMatrix object
% sorted approximately from high-level to low-level

methods
    function matrix = getMatrix(this)
        % Returns the content of this block-matrix as a matlab array
        %
 
        % allocate size for result matrix
        siz = size(this);
        matrix = zeros(siz);
         
        % determine block dimensions along each dimension
        parts1 = blockPartition(this.dims, 1);
        parts2 = blockPartition(this.dims, 2);
        
        % iterate over diagonal blocks
        for iBlock = 1:min(length(parts1), length(parts2))
            block = getBlock(this, iBlock, iBlock);
            rowInds = blockIndices(parts1, iBlock);
            colInds = blockIndices(parts2, iBlock);
            matrix(rowInds, colInds) = block;
        end
    end
 
    function block = getBlock(this, row, col)
        % return the (i-th, j-th) block 
        %
        %   BLK = getBlock(BM, ROW, COL)
        %
        
        if row == col
            % extract data element corresponding to block.
            block = this.diags{row};
            
        else
            % determine row indices of block rows
            parts1 = blockPartition(this.dims, 1);
            parts2 = blockPartition(this.dims, 2);
            
            % returns a zeros matrix of the appropriate size
            block = zeros(parts1(row), parts2(col));
        end
    end
    
    function setBlock(this, row, col, blockData)
        % set the data for the (i-th, j-th) block 
        %
        %   setBlock(BM, ROW, COL, DATA)
        %   ROW and COL indices should be equal
        %
        
        % check ROW and COL equality
        if row ~= col
            error('row and column indices should be the same');
        end
        
        % determine row indices of block rows
        parts1 = getBlockDimensions(this.dims, 1);
        rowInds = blockIndices(parts1, row)';
        
        % check number of rows of input data
        if ~isscalar(blockData) && length(rowInds) ~= size(blockData, 1)
            error('block data should have %d rows, not %d', ...
                length(rowInds), size(blockData, 1));
        end

        % determine column indices of block columns
        parts2 = getBlockDimensions(this.dims, 2);
        colInds = blockIndices(parts2, col);
        
        % check number of columns of input data
        if ~isscalar(blockData) && length(colInds) ~= size(blockData, 2)
            error('block data should have %d columns, not %d', ...
                length(colInds), size(blockData, 2));
        end

        % extract data element corresponding to block. 
        this.diags{row} = blockData;
    end
    
    function block = diagonalBlock(this, ind)
        % Return the block at a specified index on the diagonal
        block = this.diags{ind};
    end
    
    function blocks = diagonalBlocks(this)
        % Return the list of blocks on the diagonal
        blocks = this.diags;
    end
end


%% Methods that depends uniquely on BlockDimensions object
methods
    function dims = blockDimensions(this, varargin)
        % Return the dimensions of the block in the specified dimension
        %
        %   DIMS = blockDimensions(BM)
        %   Returns the block-dimension of this block matrix, as a
        %   BlockDimensions object.
        %   
        %   DIMS = blockDimensions(BM, IND)
        %   Returns the BlockDimension for the specified dimension,as an
        %   instance of IntegerPartition. 
        %
        if nargin == 1
            dims = this.dims;
        else
            dim = varargin{1};
            dims = getBlockDimensions(this.dims, dim);
        end
    end
    
    function dim = dimensionality(this)
        % Return the number of dimensions of this block matrix (usually 2)
        dim = dimensionality(this.dims);
    end
    
    function siz = blockSize(this, varargin)
        % Return the number of blocks of this BlockMatrix
        %
        % N = blockSize(BM);
        % N = blockSize(BM, DIM);
        %
        siz = blockSize(this.dims, varargin{:});
    end

    function n = blockNumber(this)
        % Return the total number of blocks of this BlockMatrix
        n = prod(blockSize(this));
    end
end


%% Apply functions on inner data
methods
    function res = blockNorm(this, varargin)
        % Computes the Block-norm of this BlockDiagonal
        %
        % NORM = blockNorm(BM)
        % returns the norm as a block diagonal: the resulting block matrix
        % is a scalar block diagonal matrix (all blocks have 1 row and 1
        % column), with the same block-size as the original matrix.
        % 
        
        % compute size of result 
        nDiags = length(this.diags);
        resNorm = cell(1, nDiags);
        
        % iterate over diagonal blocks
        for i = 1:nDiags
            % compute norm of current block
            resNorm{i} = norm(this.diags{i}, varargin{:});
        end
        
        % convert to block-diagonal instance
        res = BlockDiagonal(resNorm);
    end
        
    function res = fapply(fun, this, varargin)
        % Apply any function to the inner block matrix data

        newData = cell(1, length(this.diags));
        for i = 1:length(this.diags)
            newData{i} = fun(this.diags{i}, varargin{:});
        end
        res = BlockDiagonal(newData);
    end
    
    function reveal(this)
        % Reveal the structure of the block-Diagonal Matrix in a condensed way
        
        % extract block partitions in each direction
        parts1 = blockDimensions(this, 1);
        parts2 = blockDimensions(this, 2);
        
        % display first line with size of column blocks
        pattern = ['   ' repmat('%3d', 1, length(parts2)) '\n'];
        fprintf(pattern, parts2.terms);
        
        % display each block-row, with '0' for non-diagonal blocks
        pattern = '%3d%s\n';
        for iRow = 1:length(parts1)
            str1 = repmat('  0', 1, max(0, iRow - 1));
            str2 = repmat('  0', 1, max(0, length(parts1) - iRow));
            fprintf(pattern, parts1(iRow), [str1 '  +' str2]);
        end
    end

end


%% Overload some native methods
methods
    function siz = size(this, varargin)
        % Return the size in each direction of this block matrix object
        siz = size(this.dims, varargin{:});
    end
    
    function res = transpose(this)
        % transpose this BlockDiagonal Matrix
        res = ctranspose(this);
    end
    
    function res = ctranspose(this)
        % overload the transpose operator for BlockDiagonal object
        
        % Transpose each block
        nDiags = length(this.diags);
        diags2 = cell(nDiags, 1);
        for i = 1:nDiags
            diags2{i} = this.diags{i}';
        end
        
        % Creates the new BlockDiagonal object (Block dimensions are
        % computed automatically in constructor)
        res = BlockDiagonal(diags2);
    end
    
    function varargout = subsasgn(this, subs, value)
        % Override subsasgn function for BlockMatrix objects
        
        % extract current indexing info
        s1 = subs(1);
        type = s1.type;
        
        if strcmp(type, '.')
            % in case of dot reference, use builtin
            
            % if some output arguments are asked, use specific processing
            if nargout > 0
                varargout = cell(1);
                varargout{1} = builtin('subsasgn', this, subs, value);
            else
                builtin('subsasgn', this, subs, value);
            end
            
        elseif strcmp(type, '()')
            % In case of parens reference, index the inner data
            error('BlockDiagonal:subsasgn', 'Can not manage parens reference');
            
        elseif strcmp(type, '{}')
            % In case of braces indexing, use blocks
            
            ns = length(s1.subs);
            if ns == 1
                % linear indexing of block
                blockRow = s1.subs{1};
                setBlock(this, blockRow, blockRow, value);
                
            elseif ns == 2
                % two indices: row and col indices of blocks should be the same
                blockRow = s1.subs{1};
                blockCol = s1.subs{2};
                if any(blockRow ~= blockCol)
                    error('row indices should match column indices');
                end
                setBlock(this, blockRow, blockRow, value);
                
            else
                error('too many indices for identifying diagonal block');
            end

            
        else
            error('BlockDiagonal:subsasgn', 'Can not manage such reference');
        end
        
        if nargout > 0
            varargout{1} = this;
        end

    end

end

end % end classdef

