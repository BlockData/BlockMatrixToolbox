classdef BlockMatrix < AbstractBlockMatrix
%BLOCKMATRIX Matrix that can be divided into several blocks
%
%   BlockMatrix objects can be constructed in several ways:
%   data = reshape(1:28, [4 7]);
%   % construction from a cell array of integer partitions
%   BM = BlockMatrix.create(data, {[2 2], [2 3 2]});
%   % construction from integer partitions in each direction
%   BM = BlockMatrix.create(data, [2 2], [2 3 2]);
%   % construction from a BlockDimension object
%   DIMS = BlockDimensions({[2 2], [2 3 2]});
%   BM = BlockMatrix.create(data, DIMS);
%
%   Example
%     data = reshape(1:28, [7 4])';
%     dims = BlockDimensions({[2 2], [2 3 2]});
%     BM = BlockMatrix.create(data, dims);
%     disp(BM);
%
%   See also
%     BlockDiagonal, BlockDimensions, AbstractBlockMatrix
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-02-19,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % contains the dimensions of blocks in each dimension, as a
    % BlockDimensions object
    dims;
    
    % contains the data, as an array. 
    data;
    
end % end properties

%% Static methods
% A collection of static methods, mostly static factories.
methods (Static)
    function BM = create(data, varargin)
        %CREATE Create a new Block Matrix
        %
        % Example
        %   BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        %   reveal(BM)
        %        2  3  2
        %     2  +  +  +
        %     2  +  +  +
        %
        BM = BlockMatrix(data, varargin{:});
    end
    
    function res = zeros(blockDims)
        %ZEROS Create an empty BlockMatrix with specified block-dimensions
        %
        % Example
        %   BD = BlockDimensions({[2 2], [2 3 2]});
        %   BM = BlockMatrix.zeros(BD);
        %   reveal(BM)
        %        2  3  2
        %     2  +  +  +
        %     2  +  +  +
        %
        
        if ~isa(blockDims, 'BlockDimensions')
            error('Requires an instance of BlockDimensions as input');
        end
        
        % create empty array with appropriate size
        array = zeros(size(blockDims));
        
        % encapsulate into BlockMatrix object
        res = BlockMatrix(array, blockDims);
    end
    
    function BM = oneBlock(mat)
        %ONEBLOCK Convert a matrix or a block-matrix to a 1-1 BlockMatrix
        %
        %   BM = BlockMatrix.oneBlock(MAT)
        %   MAT is either a standard matlab array, or a BlockMatrix object.
        %   The function converts the input MAT into a BlockMatrix object
        %   such that the result has same data, but is divided into a
        %   single block in each dimension.
        %
        %   Example
        %   BM = BlockMatrix.oneBlock(magic(3));
        %   reveal(BM)
        %         3
        %      3  +
        %
        %   See also
        %     BlockMatrix, scalarBlock, uniformBlocks
        
        
        % eventually converts to single matrix
        if isa(mat, 'AbstractBlockMatrix')
            mat = getMatrix(mat);
        end
        
        % extract dimensions
        n = size(mat, 1);
        p = size(mat, 2);

        % create new BlockMatrix object
        BM = BlockMatrix(mat, n, p);
    end
    
    function BM = scalarBlock(mat)
        %SCALARBLOCK Convert a matrix to a BlockMatrix with only scalar blocks
        %
        %   BM = BlockMatrix.scalarBlock(MAT)
        %   MAT is either a standard matlab array, or a BlockMatrix object.
        %   The function converts the input MAT into a BlockMatrix object
        %   such that the result has same data, but is divided into 1-by-1
        %   blocks in each dimension.
        %
        %   Example
        %   BM = BlockMatrix.scalarBlock(magic(3));
        %   reveal(BM)
        %          1  1  1
        %       1  +  +  +
        %       1  +  +  +
        %       1  +  +  +
        %
        %   See also
        %     BlockMatrix, oneBlock, uniformBlocks
        
        % eventually converts to single matrix
        if isa(mat, 'AbstractBlockMatrix')
            mat = getMatrix(mat);
        end
        
        % extract dimensions
        n = size(mat, 1);
        p = size(mat, 2);

        % create new BlockMatrix object
        BM = BlockMatrix(mat, ones(1, n), ones(1, p));
    end
    
    function BM = uniformBlocks(mat, blockSize)
        % Convert a (block-)matrix to a block matrix with uniform blocks 
        %
        %   BM = BlockMatrix.uniformBlocks(MAT, BLOCKSIZE)
        %   MAT is either a standard matlab array, or a BlockMatrix object.
        %   The function converts the input MAT into a BlockMatrix object
        %   such that the result has same data, but is divided into blocks
        %   that all have the same size, specified by BLOCKSIZE argument.
        %
        %   Example
        %   BM = BlockMatrix.uniformBlocks(magic(6), [2 3]);
        %   reveal(BM)
        %          3  3
        %       2  +  +
        %       2  +  +
        %       2  +  +
        %
        %   See also
        %     BlockMatrix, oneBlock, scalarBlock
        
        % eventually converts to single matrix
        if isa(mat, 'AbstractBlockMatrix')
            mat = getMatrix(mat);
        end
        
        % extract dimensions
        n = size(mat, 1);
        p = size(mat, 2);
        
        % size of blocks in each dimension
        bs1 = blockSize(1);
        bs2 = blockSize(2);
        
        % number of blocks in each dimension
        nb1 = n / bs1;
        nb2 = p / bs2;
        
        % check numbers are integer
        if nb1 ~= floor(nb1)
            error('%d should be divisible by %d', n, bs1);
        end
        if nb2 ~= floor(nb2)
            error('%d should be divisible by %d', p, bs2);
        end
        
        % create new BlockMatrix object
        BM = BlockMatrix(mat, repmat(bs1, 1, nb1), repmat(bs2, 1, nb2));
        
    end
end

%% Constructor
methods
    function this = BlockMatrix(varargin)
        % Constructor for BlockMatrix class
        %
        %   BM = BlockMatrix(DATA, DIMS);
        %   Creates a new BlockMatrix instance from a data array and a
        %   BlockDimensions object.
        %
        %   BM = BlockMatrix(DATA, DIMS1, DIMS2);
        %   Creates a new BlockMatrix instance from a data array and the
        %   block dimensions for rows and columns.
        %
        %   BM = BlockMatrix(BM0);
        %   Creates a new BlockMatrix instance from an existing BlockMatrix
        %   that can also be a BlockDiagonal object.
        %   
        %   Examples:
        %     % construction from a cell array of integer partitions
        %     data = reshape(1:28, [7 4])';
        %     BM = BlockMatrix(data, {[2 2], [2 3 2]});
        %     % construction from integer partitions in each direction
        %     BM = BlockMatrix(data, [2 2], [2 3 2]);
        %     % construction from a BlockDimension object
        %     DIMS = BlockDimensions({[2 2], [2 3 2]});
        %     BM = BlockMatrix(data, DIMS);
        %
        %
        
        % (the different cases are sorted by order of expected frequency)
        if nargin == 2
            if isnumeric(varargin{1})
                % initialisation constructor
                this.data = varargin{1};
                
            elseif isa(varargin{1}, 'AbstractBlockMatrix')
                % copy constructor
                bm = varargin{1};
                this.data = getMatrix(bm);
            else
                error('first argument must be a matrix or a block matrix');
            end
            
            % Second argument represents block dimension. It can be either
            % a BlockDimensions object, or a cell array containing the
            % block sizes in each dimension. 
            var2 = varargin{2};
            if isa(var2, 'BlockDimensions')
                this.dims = var2;
            elseif iscell(var2)
                this.dims = BlockDimensions(var2);
            else
                error('second argument must be a cell array or a BlockDimensions object');
            end

            % also checks data and block dimensions match together
            checkDimensionsValidity();
            
        elseif nargin == 3
            % get data array and check validity
            if ~isnumeric(varargin{1}) || ndims(varargin{1})~=2 %#ok<ISMAT>
                error('First argument must be a 2D numeric array');
            end
            this.data = varargin{1};
            
            % create block-dimensions from the two last arguments
            rowdims = IntegerPartition(varargin{2});
            coldims = IntegerPartition(varargin{3});
            this.dims = BlockDimensions({rowdims, coldims});
            
            % also checks data and block dimensions match together
            checkDimensionsValidity();
            
        elseif nargin == 1
            % copy constructor, from another BlockMatrix object
            if isa(varargin{1}, 'AbstractBlockMatrix')
                bm = varargin{1};
                this.data = getMatrix(bm);
                this.dims = blockDimensions(bm);
            else
                error('copy constructor requires another BlockMatrix object');
            end
            
        elseif isempty(varargin)
            % empty constructor: populate with default data
            this.data = reshape(1:28, [7 4])';
            this.dims = BlockDimensions({[2 2], [2 3 2]});
            
        else
            error('Requires two or three input arguments');
        end

        
        function checkDimensionsValidity()
            % Check that data and block dimensions match together
            
            % string pattern for error message
            pattern = 'Input data have %1$d %3$s, but block dimensions specifies %2$d %3$s';
            
            % check rows
            siz1 = size(this.data, 1);
            bdim1 = sum(this.dims{1});
            if siz1 ~= bdim1
                error(pattern, siz1, bdim1, 'rows');
            end
            
            % check columns
            siz2 = size(this.data, 2);
            bdim2 = sum(this.dims{2});
            if siz2 ~= bdim2
                error(pattern, siz2, bdim2, 'columns');
            end
        end
    end

end % end constructors


%% Methods specific to BlockMatrix object
% sorted approximately from high-level to low-level

methods
    function matrix = getMatrix(this)
        % Return the content of this block-matrix as a matlab array
        %
        % For a BlockMatrix object BM, this is equivalent to 
        % matrix = BM.data;
        %
        matrix = this.data;
    end
    
    function block = getBlock(this, row, col)
        % Return the content of the (i-th, j-th) block 
        %
        %   BLK = getBlock(BM, ROW, COL)
        %
        
        % determine row indices of block rows
        parts1 = getBlockDimensions(this.dims, 1);
        rowInds = blockIndices(parts1, row)';

        % determine column indices of block columns
        parts2 = getBlockDimensions(this.dims, 2);
        colInds = blockIndices(parts2, col);
        
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
        % Set the content of the (i-th, j-th) block to specified matrix
        %
        %   setBlock(BM, ROW_IND, COL_IND, DATA)
        %
        
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
    function dims = blockDimensions(this, varargin)
        % Return the block-dimensions of this block-matrix
        %
        %   DIMS = blockDimensions(BM)
        %   Returns the block-dimension of this block matrix, as a
        %   BlockDimensions object.
        %   
        %   DIMS = blockDimensions(BM, IND)
        %   Returns the Block Dimension for the specified dimension, as an
        %   instance of IntegerPartition. 
        %
        
        if nargin == 1
            dims = this.dims;
        else
            dim = varargin{1};
            dims = getBlockDimensions(this.dims, dim);
        end
    end

    function dims = getBlockDimensions(this, varargin)
        % deprecated: use size instead
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getBlockDimensions'' is obsolete, use ''blockDimensions'' instead');
        dims = blockDimensions(this, varargin{:});
    end
    
    function dim = dimensionality(this)
        % Return the number of dimensions of this block matrix (usually 2)
        dim = dimensionality(this.dims);
    end
    
    function siz = getSize(this, varargin)
        % Return the size in each direction of this block matrix object
        % deprecated: use size instead
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getSize'' is obsolete, use ''size'' instead');
        siz = size(this.dims, varargin{:});
    end
    
    function varargout = blockSize(this, varargin)
        % Return the number of blocks in each direction
        %
        % N = blockSize(BM);
        % N = blockSize(BM, DIM);
        % [N1, N2] = blockSize(BM);
        %
        
        if nargout <= 1
            varargout = {blockSize(this.dims, varargin{:})};
        else
            varargout = {blockSize(this.dims, 1), blockSize(this.dims, 2)};
        end
    end

    function n = blockNumber(this)
        % Return the total number of blocks of this BlockMatrix
        n = prod(blockSize(this));
    end
    
    function n = getBlockNumber(this, varargin)
        % Return the total number of blocks in this block matrix, or the
        % number of blocks in a given dimension
        %
        % deprecated: use blockSize instead
        
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getBlockNumber'' is obsolete, use ''blockSize'' instead');
        n = getBlockNumber(this.dims, varargin{:});
    end
    
    function n = getBlockNumbers(this)
        % Return the number of blocks in each dimension
        % deprecated: use blockSize instead
        
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getBlockNumbers'' is obsolete, use ''blockSize'' instead');
        n = getBlockNumbers(this.dims);
    end
end


%% Overload some arithmetic operators
methods
    function res = blockNorm(this, varargin)
        % Computes the Block-norm of this BlockMatrix
        %
        % NORM = blockNorm(BM)
        % returns the norm as a block matrix: the resulting block matrix is
        % a scalar block matrix (all blocks have 1 row and 1 column), with
        % the same block-size as the original matrix.
        % 
        
        % compute size of result (corresponding to the "block-size")
        siz = blockSize(this);
        res = BlockMatrix.scalarBlock(zeros(siz));
        
        % iterate over blocks
        for i = 1:siz(1)
            for j = 1:siz(2)
                % compute norm of current block
                blockNorm = norm(getBlock(this, i, j), varargin{:});
                setBlock(res, i, j, blockNorm);
            end
        end
    end
    
    function res = fapply(fun, this, varargin)
        % Apply any function to the inner block matrix data
        %
        % RES = fapply(FUN, BM)
        %
        % Example
        % data = reshape((1:28)', [7 4])';
        % BM = BlockMatrix(data, [2 2], [2 3 2]);
        % meanBM = fapply(@mean, BM);
        % reveal(meanBM)
        %      2  3  2
        %   1  +  +  +
        %   1  +  +  +
        %
        % meanBM2 = fapply(@(x)mean(x, 2), BM);
        % reveal(meanBM2)
        %      1  1  1
        %   2  +  +  +
        %   2  +  +  +        
        %
        
        % allocate memory for intermediate results
        blockRes = cell(blockSize(this));
        blockDims1 = zeros(1, blockSize(this, 1));
        blockDims2 = zeros(1, blockSize(this, 2));
        
        % compute result for each block, and resulting dimensions
        for iBlock = 1:blockSize(this, 1)
            for jBlock = 1:blockSize(this, 2)
                block = getBlock(this, iBlock, jBlock);
                res0 = fun(block);
                
                blockRes{iBlock, jBlock} = res0;
                blockDims1(iBlock) = size(res0, 1);
                blockDims2(jBlock) = size(res0, 2);
            end
        end
        
        % create the result block matrix
        BD = BlockDimensions({blockDims1, blockDims2});
        res = BlockMatrix.zeros(BD);
        
        % populate result with individual block results
        for iBlock = 1:blockSize(this, 1)
            for jBlock = 1:blockSize(this, 2)
                setBlock(res, iBlock, jBlock, blockRes{iBlock, jBlock});
            end
        end
        
% old version
%         newData = fun(this.data, varargin{:});
%         res = BlockMatrix(newData, this.dims);
    end
end


%% Overload array manipulation methods
methods   
    function varargout = size(this, varargin)
        % Return the size in each direction of this block matrix object
        % 
        % SIZ = size(BM);
        % SIZI = size(BM, DIR);
        % [S1, S2] = size(BM);
        % Results are equivalent as using size() on the result of
        % getMatrix(A).
        %
        % Example
        %   BM = BlockMatrix(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        %   size(BM)
        %   ans =
        %        4   7
        %
        % See also
        %   blockSize
        
        
        if nargout <= 1
            varargout = {size(this.dims, varargin{:})};
        else
            varargout = {size(this.dims, 1), size(this.dims, 2)};
        end
    end

    function res = transpose(this)
        % transpose this BlockMatrix
        res = ctranspose(this);
    end
    
    function res = ctranspose(this)
        % overload the transpose operator for BlockMatrix object
        
        % ensure the new matrix is a rectangular array, with the new size
        siz = size(this);
        data2 = reshape(this.data, siz)';
        
        % transpose the BlockDimensions object, and create new BlockMatrix
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
        data2 = reshape(this.data, size(this));
        dims2 = this.dims;
        
        for i = 1:length(varargin)
            var = varargin{i};
            
            dataToAdd = reshape(var.data, size(var));
            if size(dataToAdd, 1) ~= size(data2, 1)
                error('BlockMatrices should have same number of rows');
            end
            
            data2 = [data2 dataToAdd]; %#ok<AGROW>
            dims2 = [dims2 var.dims]; %#ok<AGROW>
        end
        
        res = BlockMatrix(data2, dims2);
    end
    
    function res = vertcat(this, varargin)
        % Override the vertical concatenation operator
        
        % initialize block dimension and data to that of first BlockMatrix
        data2 = reshape(this.data, size(this));
        dims2 = this.dims;
        
        for i = 1:length(varargin)
            var = varargin{i};
            
            dataToAdd = reshape(var.data, size(var));
            if size(dataToAdd, 2) ~= size(data2, 2)
                error('BlockMatrices should have same number of columns');
            end
            
            data2 = [data2 ; dataToAdd]; %#ok<AGROW>
            dims2 = [dims2 ; var.dims]; %#ok<AGROW>
        end
        
        res = BlockMatrix(data2, dims2);
    end
end    
    

%% Overload array indexing methods

methods
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
            
            % different processing if 1 or 2 indices are used
            ns = length(s1.subs);
            if ns == 1
                % one index: use linearised indices
                
                % check that indices are within image bound
                this.data(s1.subs{:});
                
                this.data(s1.subs{1}) = value;
                
            elseif ns == 2
                % two indices: parse row and column indices
                
                % check that indices are within image bound
                this.data(s1.subs{:});
                
                % extract corresponding data
                this.data(s1.subs{:}) = value;
           
            else
                error('BlockMatrix:subsasgn', 'Too many indices');
            end
            
        elseif strcmp(type, '{}')
            % In case of braces indexing, use blocks
            
            ns = length(s1.subs);
            if ns == 2
                % Computes row and col indices of blocks
                rowBlockInds = s1.subs{1};
                if ischar(rowBlockInds) && strcmp(rowBlockInds, ':')
                    rowBlockInds = 1:blockSize(this, 1);
                end
                colBlockInds = s1.subs{2};
                if ischar(colBlockInds) && strcmp(colBlockInds, ':')
                    colBlockInds = 1:blockSize(this, 2);
                end
                
                
                % extract block partitions in each direction
                parts1 = blockDimensions(this, 1);
                parts2 = blockDimensions(this, 2);
                rowInds = blockIndices(parts1, rowBlockInds);
                colInds = blockIndices(parts2, colBlockInds);

                % process the case of chained subsref
                if length(subs) > 1
                    s2 = subs(2);
                    ns2 = length(s2.subs);
                    
                    % can only handle parens indexing
                    if ~strcmp(subs(2).type, '()')
                        error('Requires parens indexing as second reference');
                    end
                    
                    % manages two subs for the moment
                    if ns2 == 2
                        rowInds = rowInds(s2.subs{1});
                        colInds = colInds(s2.subs{2});
                    else
                        error('Requires two indices for identifying subblocks');
                    end
                end
                
                % fill sub matrix with input values
                this.data(rowInds, colInds) = value;
            else
                error('Requires two indices for identifying blocks');
            end

            
        else
            error('BlockMatrix:subsasgn', 'Can not manage such reference');
        end
        
        if nargout > 0
            varargout{1} = this;
        end

    end
end

end % end classdef

