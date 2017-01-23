classdef BlockMatrix < BlockTensor
%BLOCKMATRIX  Abstract parent class of DenseBlockMatrix and BlockDiagonal 
%
%   Class BlockMatrix
%
%   Example
%     data = reshape(1:28, [7 4])';
%     BM = BlockMatrix.create(data, [1 2 1], [2 3 2]);
%
%   See also
%     DenseBlockMatrix, BlockDiagonal
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-04-22,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


% no properties are defined, nor constructor

%% Interface Methods
% Methods in this cell are declared for implementation in sub-classes

methods (Abstract)
    % Return the content of this block-matrix as a matlab array
    matrix = getMatrix(this)
    
    % Return the dimensions of the block in the specified dimension
    dims = blockDimensions(this, dim)
    
    % Return the number of dimensions of this block matrix (usually 2)
    dim = dimensionality(this)
       
    % Return the total number of blocks in this block matrix
    n = blockNumber(this, varargin)
    
    % Return the number of blocks in each dimension
    n = getBlockNumbers(this)
    
    % Return the number of blocks in each dimension
    n = blockSize(this, varargin)
    
    % Returns the block content at a given position
    block = getBlock(this, row, col)
    
    % Updates the block content at a given position
    block = setBlock(this, row, col, blockContent)

end % end abstract methods

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
        BM = DenseBlockMatrix(data, varargin{:});
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
        
        % encapsulate into DenseBlockMatrix object
        res = DenseBlockMatrix(array, blockDims);
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
        if isa(mat, 'BlockMatrix')
            mat = getMatrix(mat);
        end
        
        % extract dimensions
        n = size(mat, 1);
        p = size(mat, 2);

        % create new DenseBlockMatrix object
        BM = DenseBlockMatrix(mat, n, p);
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
        if isa(mat, 'BlockMatrix')
            mat = getMatrix(mat);
        end
        
        % extract dimensions
        n = size(mat, 1);
        p = size(mat, 2);

        % create new DenseBlockMatrix object
        BM = DenseBlockMatrix(mat, ones(1, n), ones(1, p));
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
        if isa(mat, 'BlockMatrix')
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
        BM = BlockMatrix.create(mat, repmat(bs1, 1, nb1), repmat(bs2, 1, nb2));
        
    end
end


%% Computation Methods
% Some methods are performed that only rely on abstract methods.

methods
    function res = blockProduct(this, that, type)
        % Compute newly defined block-matrix product of two block-matrices
        %
        % RES = blockProduct(BM1, BM1, TYPE)
        % BM1 and BM2 are two block matrices, and TYPE is a string
        % representing the type of blck product, for eaxmple 'uu', 'hh',
        % 'ks'...
        %
        % Some conditions exist on the block-dimensions of input block
        % matrices, depending on the type of block-product applied.
        %
        % Example
        %   BM = BlockMatrix.create(reshape(1:28, [4 7]), [1 2 1], [2 3 2]);
        %   res = blockProduct(BM, BM', 'uu');
        %   disp(res)
        %   DenseBlockMatrix object with 4 rows and 4 columns
        %     row dims: 2 2
        %     col dims: 2 2
        % 
        %    140.0000    336.0000    532.0000    728.0000   
        %    336.0000    875.0000    1.41e+03    1.95e+03   
        %    532.0000    1.41e+03    2.30e+03    3.18e+03   
        %    728.0000    1.95e+03    3.18e+03    4.40e+03   

        
        % check validity of TYPE argument
        if ~ischar(type)
            error('the TYPE argument must be a char array');
        end
        validTypes = {...
            'uu', 'uh', 'uk', 'us', ...
            'hu', 'hh', 'hk', 'hs', ...
            'ku', 'kh', 'kk', 'ks', ...
            'su', 'sh', 'sk', 'ss'};
        if ~ismember(type, validTypes)
            error('the TYPE argument is not a valid string');
        end
        
        switch lower(type(1))

            case 's'
                % scalar product along blocks
                % we need one scalar
                switch lower(type(2))
                    case 's'
                        res = blockProduct_ss(this, that);
                    case 'h'
                        res = blockProduct_sh(this, that);
                    case 'u'
                        res = blockProduct_su(this, that);
                    case 'k'
                        res = blockProduct_sk(this, that);
                end
                
            case 'h'
                % hadamard product along blocks
                % we need: same number of blocks in each direction
                switch lower(type(2))
                    case 's'
                        res = blockProduct_hs(this, that);
                    case 'h'
                        res = blockProduct_hh(this, that);
                    case 'u'
                        res = blockProduct_hu(this, that);
                    case 'k'
                        res = blockProduct_hk(this, that);
                end
                
            case 'u'
                % usual product along blocks
                % we need: blockdims1(2) == blockdims2(1)
                switch lower(type(2))
                    case 's'
                        res = blockProduct_us(this, that);
                    case 'h'
                        res = blockProduct_uh(this, that);
                    case 'u'
                        res = blockProduct_uu(this, that);
                    case 'k'
                        res = blockProduct_uk(this, that);
                end
                
            case 'k'
                % kroenecker product along blocks
                % we need: (what ?)
                switch lower(type(2))
                    case 's'
                        res = blockProduct_sk(this, that);
                    case 'h'
                        res = blockProduct_hk(this, that);
                    case 'u'
                        res = blockProduct_uk(this, that);
                    case 'k'
                        error('The ''kk''-type product is not implemented');
                end
                
        end
    end
    
    function nr = blockRows(this)
        % Return the number of block rows
        %
        % Example:
        % BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        % blockRows(BM)
        % ans = 
        %     2
        %
        % See also
        %    blockSize, blockCols
        
        nr = blockSize(this, 1);
    end
    
    function nc = blockCols(this)
        % Return the number of block columns
        %
        % Example:
        % BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        % blockCols(BM)
        % ans = 
        %     3
        %
        % See also
        %    blockSize, blockRows
        
        nc = blockSize(this, 2);
    end
    
    function res = blockTranspose(this)
        % Return a new block-matrix with each block transposed
        %
        % Example
        %   rng(100);
        %   BM = BlockMatrix.create(randi(20, 6, 6), [3 3], [2 2 2]);
        %   blockTranspose(BM)
        %
        %     ans = 
        %     DenseBlockMatrix object with 4 rows and 9 columns
        %       row dims: 2 2
        %       col dims: 3 3 3
        % 
        %        11       6       9       4       3       5       4       8       1   
        %        14      17       3      17       6       9      12      13       3   
        %        17       1       3      20      17       4       6      16       1   
        %        12      18       5      19      17       7       8       1      18   
       
        
        % first check validity of input matrix
        BD = blockDimensions(this);
        if ~isUniform(BD{1}) || ~isUniform(BD{2})
            error('requires a block-uniform DenseBlockMatrix');
        end
        
        % compute size of result
        blockSize1 = BD{1}(1);
        blockSize2 = BD{2}(1);
        parts1 = repmat(blockSize2, 1, length(BD{1}));
        parts2 = repmat(blockSize1, 1, length(BD{2}));
        BD2 = BlockDimensions({parts1, parts2});
        
        % create result block matrix
        res = BlockMatrix.zeros(BD2);

        % transpose each block
        for iBlock = 1:blockSize(this, 1)
            for jBlock = 1:blockSize(this, 2)
                setBlock(res, iBlock, jBlock, getBlock(this, iBlock, jBlock)');
            end
        end
    end
end


%% Test block partition types

methods
    function tf = isOneBlock(this)
        % Check if a BlockMatrix is divided in one block in each direction
        tf = all(blockSize(this) == [1 1]);
    end
    
    function tf = isScalarBlock(this)
        % Check if a BlockMatrix is divided in 1-1 blocks in each direction
        tf = all(blockSize(this) == size(this));
    end
    
    function tf = isUniformBlock(this)
        % Check if a BlockMatrix is divided in blocks with all the same size
        unif1 = isUniform(blockDimensions(this, 1));
        unif2 = isUniform(blockDimensions(this, 2));
        tf = unif1 && unif2;
    end
    
    function tf = isVectorBlock(this)
        % Check if a BlockMatrix is divided in 1 block in at least one direction
        tf = any(blockSize(this) == 1);
    end
    
    function tf = isPositiveDefinite(this)
        % Test is a BlockMatrix is symmetric positive definite
        %
        % Uses Sylvester's criterion. Throws an error if matrix is not
        % square and symmetric.
        
        % first check square size
        [n1, n2] = size(this);
        if n1 ~= n2
            error('Criterion is valid only for square matrices');
        end
        
        % check symmetric matrix
        mat = getMatrix(this);
        if any(any(mat ~= mat'))
            error('Criterion is valid only for symmetric matrices');
        end
        
        % compute all determinants and check their sign
        for i = 1:n1
            sub = mat(1:i, 1:i);
            if det(sub) < 0
                tf = false;
                return;
            end
        end
        
        % if no test was escaped, then all determinants are positive, 
        % and the matrix is positive definite
        tf = true;
    end
    
    function tf = isSymmetric(this)
        % Test is a BlockMatrix is symmetric
        %
        % Throws an error if matrix is not square.
        
        % first check square size
        [n1, n2] = size(this);
        if n1 ~= n2
            error('Requires a square matrix');
        end
        
        % check symmetric matrix
        mat = getMatrix(this);
        tf = ~any(any(mat ~= mat'));
    end
end


%% Extract sub matrix
methods
    function res = subMatrix(this, rowBlockInds, colBlockInds)
        % Return a new Block-Matrix corresponding to specified block indices
        %
        % BM2 = subMatrix(BM, RBINDS, CBINDS);
        %
        % Example
        %   BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        %   BM2 = subMatrix(BM, 1, [2 3]);
        %   ans = 
        %   DenseBlockMatrix object with 2 rows and 5 columns
        %       row dims: 2
        %       col dims: 3 2
        %             3           4           5           6           7   
        %            10          11          12          13          14  
        %
        % See Also
        %   subsref, getBlock
        
        % get integer partitions for each dimension
        parts1 = blockDimensions(this, 1);
        parts2 = blockDimensions(this, 2);
        
        % compute new block dimensions
        newParts1 = parts1{rowBlockInds};
        newParts2 = parts2{colBlockInds};
        newDims = BlockDimensions({newParts1, newParts2});

        res = BlockMatrix.zeros(newDims);
        for iRow = 1:length(newParts1)
            % index of block-row in original block-matrix
            rowInd = rowBlockInds(iRow);
            for iCol = 1:length(newParts2)
                blockData = getBlock(this, rowInd, colBlockInds(iCol));
                setBlock(res, iRow, iCol, blockData);
            end
        end
    end
    
    function setSubMatrix(this, rowBlockInds, colBlockInds, newData)
        % Return a new Block-Matrix corresponding to specified block indices
        %
        % Example
        %     BM = BlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2]);
        %     BM{2, [1 3]} = ones(2, 4)
        %     BM = 
        %     DenseBlockMatrix object with 4 rows and 7 columns
        %       row dims: 2 2
        %       col dims: 2 3 2
        %             1           2           3           4           5           6           7   
        %             8           9          10          11          12          13          14   
        %             1           1          17          18          19           1           1   
        %             1           1          24          25          26           1           1   
        % 
        % see also
        %   subMatrix, subsasgn
                
        % get integer partitions for each dimension
        parts1 = blockDimensions(this, 1);
        parts2 = blockDimensions(this, 2);
        
        % compute new block dimensions
        newParts1 = parts1{rowBlockInds};
        newParts2 = parts2{colBlockInds};
        newDims = BlockDimensions({newParts1, newParts2});
        
        newBlockData = BlockMatrix.create(newData, newDims);

        for iRow = 1:length(newParts1)
            % index of block-row in original block-matrix
            rowInd = rowBlockInds(iRow);
            for iCol = 1:length(newParts2)
                blockData = getBlock(newBlockData, iRow, iCol);
                setBlock(this, rowInd, colBlockInds(iCol), blockData);
            end
        end
    end
end


%% Overload EXPONENT and LOGARITHMS methods

methods
    function res = exp(this)
        % Compute the exponential for each element of the BlockMatrix
        res = BlockMatrix.create(exp(getMatrix(this)), blockDimensions(this));
    end
    
    function res = expm1(this)
        % Compute exp(x)-1, compensating for the roundoff in exp(x) for BlockMatrix
        res = BlockMatrix.create(expm1(getMatrix(this)), blockDimensions(this));
    end
    
    function res = log(this)
        % Compute the secant in radians for each element of the BlockMatrix
        res = BlockMatrix.create(log(getMatrix(this)), blockDimensions(this));
    end
    
    function res = log10(this)
        % Compute the base 10 logarithm of each element of the BlockMatrix
        res = BlockMatrix.create(log10(getMatrix(this)), blockDimensions(this));
    end
    
    function res = logp1(this)
        % Compute log(1+x), compensating for the roundoff in 1+x.of the BlockMatrix
        res = BlockMatrix.create(logp1(getMatrix(this)), blockDimensions(this));
    end
    
    function res = log2(this)
        % Compute the base 2 logarithm of the elementsof the BlockMatrix
        res = BlockMatrix.create(log2(getMatrix(this)), blockDimensions(this));
    end
    
    function res = nextpow2(this)
        % Compute the exponents for the smallest powers of two that satisfy for each element of the BlockMatrix
        res = BlockMatrix.create(nextpow2(getMatrix(this)), blockDimensions(this));
    end
    
    function res = nthroot(this)
        % Compute the real nth root of the elements of the BlockMatrix
        res = BlockMatrix.create(nthroot(getMatrix(this)), blockDimensions(this));
    end
    
    function res = pow2(this)
        % Compute a BlockMatrix X whose elements are 2 raised to the power of each element of the BlockMatrix.
        res = BlockMatrix.create(pow2(getMatrix(this)), blockDimensions(this));
    end
    
    function res = reallog(this)
        % Compute the natural logarithm of each element of the BlockMatrix
        res = BlockMatrix.create(reallog(getMatrix(this)), blockDimensions(this));
    end
    
    function res = realsqrt(this)
        % Compute the square root of each element of the BlockMatrix
        res = BlockMatrix.create(realsqrt(getMatrix(this)), blockDimensions(this));
    end
    
    function res = sqrt(this)
        % Compute the square root of each element of the BlockMatrix
        res = BlockMatrix.create(sqrt(getMatrix(this)), blockDimensions(this));
    end
end


%% Overload some arithmetic methods
methods
    function res = sin(this)
        % Compute the sin for each element of the BlockMatrix
        res = BlockMatrix.create(sin(getMatrix(this)), blockDimensions(this));
    end

    function res = cos(this)
        % Compute the cos for each element of the BlockMatrix
        res = BlockMatrix.create(cos(getMatrix(this)), blockDimensions(this));
    end

    function res = tan(this)
        % Compute the sin for each element of the BlockMatrix
        res = BlockMatrix.create(tan(getMatrix(this)), blockDimensions(this));
    end

     function res = sind(this)
        % Compute the sind for each element of the BlockMatrix
        res = BlockMatrix.create(sind(getMatrix(this)), blockDimensions(this));
    end

    function res = cosd(this)
        % Compute the cosd for each element of the BlockMatrix
        res = BlockMatrix.create(cosd(getMatrix(this)), blockDimensions(this));
    end

   function res = tand(this)
        % Compute the tand for each element of the BlockMatrix
        res = BlockMatrix.create(tand(getMatrix(this)), blockDimensions(this));
    end

    function res = csch(this)
        % Compute the hyperbolic cosecant for each element of the BlockMatrix
        res = BlockMatrix.create(csch(getMatrix(this)), blockDimensions(this));
    end
    
    function res = acsch(this)
        % Compute the inverse hyperbolic cosecant for each element of the BlockMatrix
        res = BlockMatrix.create(acsch(getMatrix(this)), blockDimensions(this));
    end
    
    function res = sec(this)
        % Compute the secant in radians for each element of the BlockMatrix
        res = BlockMatrix.create(sec(getMatrix(this)), blockDimensions(this));
    end
    
    function res = secd(this)
        % Compute the secant in degrees for each element of the BlockMatrix
        res = BlockMatrix.create(secd(getMatrix(this)), blockDimensions(this));
    end
    
    function res = asec(this)
        % Compute the inverse secant in radians for each element of the BlockMatrix
        res = BlockMatrix.create(asec(getMatrix(this)), blockDimensions(this));
    end
    
    function res = asecd(this)
        % Compute the inverse secant in degrees for each element of the BlockMatrix
        res = BlockMatrix.create(asecd(getMatrix(this)), blockDimensions(this));
    end
    
    function res = sech(this)
        % Compute the returns the hyperbolic secant for each element of the BlockMatrix
        res = BlockMatrix.create(sech(getMatrix(this)), blockDimensions(this));
    end
    
    function res = asech(this)
        % Compute the inverse hyperbolic secant for each element of the BlockMatrix
        res = BlockMatrix.create(asech(getMatrix(this)), blockDimensions(this));
    end
    
    function res = cot(this)
        % Compute the returns the cotangent in radians for each element of the BlockMatrix
        res = BlockMatrix.create(cot(getMatrix(this)), blockDimensions(this));
    end
    
    function res = cotd(this)
        % Compute the returns the cotangent in degrees for each element of the BlockMatrix
        res = BlockMatrix.create(cotd(getMatrix(this)), blockDimensions(this));
    end
    
    function res = acot(this)
        % Compute the returns the inverse cotangent in radius for each element of the BlockMatrix
        res = BlockMatrix.create(acot(getMatrix(this)), blockDimensions(this));
    end
    
    function res = acotd(this)
        % Compute the returns the inverse cotangent in degrees for each element of the BlockMatrix
        res = BlockMatrix.create(acotd(getMatrix(this)), blockDimensions(this));
    end
    
    function res = acoth(this)
        % Compute the inverse hyperbolic cotangent for each element of the BlockMatrix
        res = BlockMatrix.create(acoth(getMatrix(this)), blockDimensions(this));
    end
    
    function res = coth(this)
        % Compute the hyperbolic cotangent  for each element of the BlockMatrix
        res = BlockMatrix.create(coth(getMatrix(this)), blockDimensions(this));
    end
end

%% Overload some arithmetic methods

methods
    function res = norm(this, varargin)
        % Compute the norm of this DenseBlockMatrix
        %
        % NORM = norm(BM)
        % returns the norm as a scalar
        % 
        res = norm(getMatrix(this), varargin{:});
    end

    function res = plus(this, that)
        % Addition of two block-matrices
        % 
        % usage:
        %   X = A + B
        
        sizA = size(this);
        sizB = size(that);
        if any(sizA ~= sizB) && ~(prod(sizA)==1 || prod(sizB)==1)
            error('Both inputs must have same size, or one must be scalar');
        end
        
        dataA = this;
        dataB = that;
        parent = [];
        
        if isa(this, 'BlockMatrix')
            dataA = this.data;
            parent = this;
        end
        
        if isa(that, 'BlockMatrix')
            dataB = that.data;
            parent = that;
        end
       
        % compute data of resulting matrix
        dataX = dataA + dataB;
        
        % create resulting object
        dimsX = blockDimensions(parent);
        res = BlockMatrix.create(dataX, dimsX);
    end
    
    function res = minus(this, that)
        % Subtraction of two block-matrices
        % 
        % usage:
        %   X = A - B
        
        sizA = size(this);
        sizB = size(that);
        if any(sizA ~= sizB) && ~(prod(sizA)==1 || prod(sizB)==1)
            error('Both inputs must have same size, or one must be scalar');
        end
        
        dataA = this;
        dataB = that;
        parent = [];
        
        if isa(this, 'BlockMatrix')
            dataA = getMatrix(this);
            parent = this;
        end
        
        if isa(that, 'BlockMatrix')
            dataB = getMatrix(that);
            if isempty(parent)
                parent = that;
            end
        end
       
        % compute data of resulting matrix
        dataX = dataA - dataB;
        
        % create resulting object
        dimsX = blockDimensions(parent);
        res = BlockMatrix.create(dataX, dimsX);
    end
    
    function res = times(this, that)
        % Element-wise multiplication of two block-matrices
        % 
        % usage:
        %   X = A .* B
        
        sizA = size(this);
        sizB = size(that);
        if any(sizA ~= sizB) && ~(prod(sizA)==1 || prod(sizB)==1)
            error('Both inputs must have same size, or one must be scalar');
        end
        
        dataA = this;
        dataB = that;
        parent = [];
        
        if isa(this, 'BlockMatrix')
            dataA = this.data;
            parent = this;
        end
        
        if isa(that, 'BlockMatrix')
            dataB = that.data;
            parent = that;
        end
       
        % compute data of resulting matrix
        dataX = dataA .* dataB;
        
        % create resulting object
        dimsX = blockDimensions(parent);
        res = BlockMatrix.create(dataX, dimsX);
    end
    
    function res = rdivide(this, that)
        % Division of two block-matrices
        % 
        % usage:
        %   X = A ./ B
        
        sizA = size(this);
        sizB = size(that);
        if any(sizA ~= sizB) && ~(prod(sizA)==1 || prod(sizB)==1)
            error('Both inputs must have same size, or one must be scalar');
        end
        
        dataA = this;
        dataB = that;
        parent = [];
        
        if isa(this, 'BlockMatrix')
            dataA = this.data;
            parent = this;
        end
        
        if isa(that, 'BlockMatrix')
            dataB = that.data;
            parent = that;
        end
       
        % compute data of resulting matrix
        dataX = dataA ./ dataB;
        
        % create resulting object
        dimsX = blockDimensions(parent);
        res = BlockMatrix.create(dataX, dimsX);
    end
    
    function res = mtimes(this, that)
        % Multiply two instances of DenseBlockMatrix
        % 
        % usage:
        %   X = A * B
        
        % check case of multiplication by a numeric value
        if isnumeric(this)
            if isscalar(this)
                res = times(this, that);
                return;
            else
                % convert the numeric array to one-block matrix to perform
                % block-matrix product
                this = BlockMatrix.oneBlock(this);
            end
        end
        if isnumeric(that)
            if isscalar(that)
                res = times(this, that);
                return;
            else
                % convert the numeric array to one-block matrix to perform
                % block-matrix product
                that = BlockMatrix.oneBlock(that);
            end
        end
        
        % get block dimensions of each matrix
        dimsA = this.dims;
        dimsB = that.dims;
        
        % check dimensionality
        if dimensionality(dimsB) ~= 2
            error('Requires another BlockTensor of dimensionality 2');
        end
        
        % total number of elements should match
        if size(dimsA, 2) ~= size(dimsB, 1)
            error('number of columns of first matrix (%d) should match number of rows of second matrix (%d)', ...
                size(dimsA, 2), size(dimsB, 1));
        end
        if blockSize(dimsA, 2) ~= blockSize(dimsB, 1)
            error('number of block columns of first matrix (%d) should match number of block rows of second matrix (%d)', ...
                blockSize(dimsA, 2), blockSize(dimsB, 1));
        end

        % compute block dimension of the resulting block-matrix
        dimsC = BlockDimensions([dimsA.parts(1) dimsB.parts(2)]);
        
        % allocate memory for result
        res = BlockMatrix.zeros(dimsC);

        % number of blocks to iterate
        nBlocks = blockSize(dimsA, 2);

        for iRow = 1:blockSize(dimsA, 1)
            for iCol = 1:blockSize(dimsB, 2)
                % Compute block (iRow, iCol), by iterating over i-th row of
                % first matrix, and j-th column of second matrix
                
                % allocate memory for current result block
                block = zeros([dimsA{1}(iRow) dimsB{2}(iCol)]);
                
                % iterate over columns of first matrix, and rows of second
                % matrix
                for iBlock = 1:nBlocks
                    blockA = getBlock(this, iRow, iBlock);
                    blockB = getBlock(that, iBlock, iCol);
                    block = block + blockA * blockB;
                end
                
                % store current block in result block matrix
                setBlock(res, iRow, iCol, block);
            end
        end
    end

end % end methods


%% Array concatenation
methods
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
        data2 = reshape(getMatrix(this), size(this));
        dims2 = blockDimensions(this);
        
        for i = 1:length(varargin)
            var = varargin{i};
            
            dataToAdd = reshape(getMatrix(var), size(var));
            if size(dataToAdd, 1) ~= size(data2, 1)
                error('BlockMatrices should have same number of rows');
            end
            
            data2 = [data2 dataToAdd]; %#ok<AGROW>
            dims2 = [dims2 var.dims]; %#ok<AGROW>
        end
        
        res = BlockMatrix.create(data2, dims2);
    end
    
    function res = vertcat(this, varargin)
        % Override the vertical concatenation operator
        
        % initialize block dimension and data to that of first BlockMatrix
        data2 = reshape(getMatrix(this), size(this));
        dims2 = blockDimensions(this);
        
        for i = 1:length(varargin)
            var = varargin{i};
            
            dataToAdd = reshape(getMatrix(var), size(var));
            if size(dataToAdd, 2) ~= size(data2, 2)
                error('BlockMatrices should have same number of columns');
            end
            
            data2 = [data2 ; dataToAdd]; %#ok<AGROW>
            dims2 = [dims2 ; var.dims]; %#ok<AGROW>
        end
        
        res = BlockMatrix.create(data2, dims2);
    end
end


%% Overload indexing methods
methods
    function varargout = subsref(this, subs)
        % Override subsref for Block Matrices objects
        %
        % BM = DenseBlockMatrix.create(reshape(1:28, [7 4])', [2 2], [2 3 2])
        % BM(2, 3)
        % ans =
        %     10
        %
        % BM{2, 3}
        % ans = 
        %    20   21
        %    27   28
        %
        
        % extract reference type
        s1 = subs(1);
        type = s1.type;
        
        % switch between reference types
        if strcmp(type, '.')
            % in case of dot reference, use builtin subsref
            
            % check if we need to return output or not
            if nargout > 0
                % if some output arguments are asked, pre-allocate result
                varargout = cell(nargout, 1);
                [varargout{:}] = builtin('subsref', this, subs);
                
            else
                % call parent function, and eventually return answer
                builtin('subsref', this, subs);
                if exist('ans', 'var')
                    varargout{1} = ans; %#ok<NOANS>
                end
            end
            
        elseif strcmp(type, '()')
            % Process parens indexing -> return elements of data matrix
            matrix = getMatrix(this);
            varargout{1} = subsref(matrix, subs);
            
        elseif strcmp(type, '{}')
            % Process braces indexing -> return new Block Matrix
            % corresponding to the specified block(s)
            ns = length(s1.subs);
            if ns ~= 2
                error('Requires two indices for identifying blocks');
            end
            
            % compute integer partition of corresponding dimension
            rowBlockInds = s1.subs{1};
            if ischar(rowBlockInds) && strcmp(rowBlockInds, ':')
                rowBlockInds = 1:blockSize(this, 1);
            end
            colBlockInds = s1.subs{2};
            if ischar(colBlockInds) && strcmp(colBlockInds, ':')
                colBlockInds = 1:blockSize(this, 2);
            end
            
            % extract the sub-blockmatrix corresponding to selected indices
            res = subMatrix(this, rowBlockInds, colBlockInds);
            varargout{1} = res;
            
            % eventually call recursively subsref
            if length(subs) > 1
                varargout{1} = subsref(res, subs(2:end));
            end
        end
    end
    
    function n = numArgumentsFromSubscript(this, subs, indexingContent)
        % Need to overload this to allow proper braces indexing
        n = numel(this);
    end
end

%% Display methods

methods
    function reveal(this)
        % Reveal the structure of the block-Matrix in a condensed way
        
        % extract block partitions in each direction
        parts1 = blockDimensions(this, 1);
        parts2 = blockDimensions(this, 2);
        
        % computes number of digits for each column
        % (use same number for all columns)
        digitsByCol = ceil(log10(parts2.terms));
        nDigits =  max(2, max(digitsByCol));
        
        colPattern = ['%' num2str(nDigits + 1) 'd'];
        pattern = ['   ' repmat(colPattern, 1, length(parts2)) '\n'];
        fprintf(pattern, parts2.terms);
        
        colPattern = [repmat(' ', 1, nDigits) '+'];
        pattern = ['%3d' repmat(colPattern, 1, length(parts2)) '\n'];
        for iRow = 1:length(parts1)
            fprintf(pattern, parts1(iRow));
        end
    end
    
    function disp(this)
        % Display the content of this DenseBlockMatrix object
        
        % loose format: display more empty lines
        isLoose = strcmp(get(0, 'FormatSpacing'), 'loose');
        
        % get DenseBlockMatrix total size
        dim = size(this);
        nRows = dim(1);
        nCols = dim(2);
        
        % Display information on block sizes
        className = class(this);
        disp(sprintf('%s object with %d rows and %d columns', className, nRows, nCols)); %#ok<DSPS>
        parts1 = blockPartition(this.dims, 1);
        disp(sprintf('  row dims: %s', formatParts(parts1))); %#ok<DSPS>
        parts2 = blockPartition(this.dims, 2);
        disp(sprintf('  col dims: %s', formatParts(parts2))); %#ok<DSPS>

        if nCols < 20 && nRows < 50
            if isLoose
                fprintf('\n');
            end
            displayData(this);
        end
        
        if isLoose
            fprintf('\n');
        end
        
        function string = formatParts(parts)
            % display parts as a list of ints, or as empty
            if isempty(parts)
                string = '(empty)';
            else
                pattern = strtrim(repmat(' %d', 1, length(parts)));
                string = sprintf(pattern, parts.terms);
            end
        end
    end
    
    
    function displayBlocks(this)
        % Display inner blocks of block matrix object
        nRowBlocks = blockSize(this.dims, 1);
        nColBlocks = blockSize(this.dims, 2);
        
        for row = 1:nRowBlocks
            for col = 1:nColBlocks
                disp(sprintf('Block (%d,%d)', row, col)); %#ok<DSPS>
                disp(getBlock(this, row, col));
            end
        end
        
    end
    
    function displayData(this)
        % Display data using different style for aternating blocks
        
        % get Block dimensions in each dimensions,
        % for the moment as a row vector of positive integers,
        % later as an instance of IntegerPartition
        dims1 = blockDimensions(this, 1);
        dims2 = blockDimensions(this, 2);
        
        % get DenseBlockMatrix total size
        nRowBlocks = length(dims1);
        nColBlocks = length(dims2);
        
        % define printing styles for alternating blocs
        styleList = {'*blue', [.7 0 0]};
        
        % the style of first block, and of 'odd' blocks
        iStyle = 1;
        
        % iterate over block-rows
        for iBlock = 1:nRowBlocks
            
            % iterate over the rows of current row of blocks
            for iRow = 1:dims1(iBlock)
                
                % style of first block
                iBlockStyle = iStyle;
                
                % iterate over blocks of current block row
                for jBlock = 1:nColBlocks
                    % extract data of current row within current block
                    blockData = getBlock(this, iBlock, jBlock);
                    rowData = blockData(iRow, :);
                    
                    % the string to display
                    stringArray = formatArray(rowData);
                    string = [strjoin(stringArray, '   ') '   '];
                    
                    % choose appropriate style
                    style = styleList{iBlockStyle};
                    
                    % display in color
                    cprintf(style, string);
                    % alternate the style of next block(-column)
                    iBlockStyle = 3 - iBlockStyle;
                end
                
                fprintf('\n');
                
            end % end iteration of rows within block-row
            
            % alternate the style for next block-row
            iStyle = 3 - iStyle;
            
        end  % end block-row iteration
        
        function stringArray = formatArray(array)
            % Convert a numerical array to a formatted cell array of strings
            stringArray = cell(size(array));
            
            for i = 1:numel(array)
                num = array(i);
                
                % choose formatting style
                if num == 0 || round(num) == num
                    fmt = '%d';
                elseif abs(num) < 1e4 && abs(num) > 1e-2
                    fmt = '%.3f';
                else
                    fmt = '%4.2e';
                end
                
                % ensure 9 digits are used, and align to the right
                stringArray{i} = sprintf('%5s', num2str(num, fmt));
            end
        end
    end % end displayData method
    
end % end display methods

end % end classdef
