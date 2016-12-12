classdef BlockDimensions < handle
%BLOCKDIMENSIONS  Store the block dimensions of a BlockMatrix data structure
%
%   A BlockDimensions object represents the partitions into blocks of a
%   block-matrix object.
%   A BlockDimensions object is represented by a series of IntegerPartition
%   object, each IntegerPartition corresponding to the block-partition
%   along a given dimension.
%
%
%   Example
%     % create a new BlockDimension object from integer partitions
%     p1 = IntegerPartition([2, 2]);
%     p2 = IntegerPartition([2, 3, 2]);
%     BD = BlockDimensions({p1, p2});
%
%     % Alternative creation syntax
%     BD = BlockDimensions({[2 2], [2, 3, 2]})
%     BD = 
%     BlockDimensions object with 2 dimensions
%     ( (2, 2), (2, 3, 2) )
%
%     % Returns the block partition along second direction
%     BD2 = BD{2}
%     BD2 = 
%     IntegerPartition object with 3 terms
%         (2, 3, 2)
%   
%     % computes the block dimension of the transposed block-matrix
%     BDT = BlockDimensions({BD{2}, BD{1}})
%     BDT = 
%     BlockDimensions object with 2 dimensions
%     ( (2, 3, 2), (2, 2) )
%
%     % returns the third block size in second direction
%     BD{2}(3)
%     ans =
%          2
%
%   See also
%     BlockMatrix, IntegerPartition

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-02-20,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % contains the dimensions of blocs in each dimension, as a cell array
    % containing instances of IntegerPartition
    parts;
    
end % end properties


%% Constructor
methods
    function this = BlockDimensions(varargin)
        % Constructor for BlockDimensions class
        %
        %   BD = BlockDimensions(PARTS);
        %   PARTS is a cell array, containing for each dimension the sizes
        %   of the blocks in this dimension.
        % 
        %   BD = BlockDimensions(PARTS, DIMS);
        %   Alternative construction of BlockDimension object. PARTS is a
        %   row vector containing all  the block dimensions, and DIMS is 
        % 
        %
        %   Example
        %     % Creates a BlockDimensions object for a BlockMatrix, divided
        %     % into two blocks in dimension 1 and into three blocks in 
        %     % dimension 2.
        %     BD = BlockDimensions({[2 2], [2 3 2]});
        %
        %     % Construction of the same BlockDimension object, using a
        %     % list of partitions and a list of term numbers
        %     BD = BlockDimensions([2 2 2 3 2], [2 3]);
        %
        
        if iscell(varargin{1})
            % Processes following cases:
            % * a cell array of IntegerPartition
            % * a cell array of integer arrays
        
            var1 = varargin{1};
            this.parts = cell(1, length(var1));
            for i = 1:length(var1)
                if isa(var1{i}, 'IntegerPartition')
                    this.parts{i} = var1{i};
                else
                    % convert integer array to IntegerPartition object
                    this.parts{i} = IntegerPartition(var1{i});
                end
            end
            
        elseif nargin == 2 && isnumeric(varargin{1}) && isnumeric(varargin{2})
            % constructor from partition array and dimension array
            allParts = varargin{1};
            dims = varargin{2};
            nd = length(dims);
            this.parts = cell(1, nd);
            ind = 1;
            for i = 1:nd
                terms = allParts(ind:(ind+dims(i)-1));
                this.parts{i} = IntegerPartition(terms);
                ind = ind + dims(i);
            end
            
        elseif isa(varargin{1}, 'BlockDimensions')
            % copy constuctor
            bd = varargin{1};
            this.parts = bd.parts;
        else
            error('input argument must be a cell array');
        end
    end

end % end constructors


%% Methods

methods
    function p = partitions(this)
        % return the vector containing all the block partitions
        % see also: dimensions
        
        % dimensionality of this BlockDimensions
        nd = length(this.parts);
        
        % compute total number of partitions
        np = 0;
        for i = 1:nd
            np = np + length(this.parts{i});
        end
        
        p = zeros(1, np);
        ind = 1;
        for i = 1:nd
            length_i = length(this.parts{i});
            p(ind:ind+length_i-1) = this.parts{i}.terms;
            ind = ind + length_i;
        end
    end
    
    function p = getPartitions(this)
        % deprecated: use "partitions" instead
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getPartitions'' is obsolete, use ''partitions'' instead');
        p = partitions(this);
    end
    
    function dims = dimensions(this)
        % return the vector containing all the dimensions
        % see also: partitions
        
        % dimensionality of this BlockDimensions
        nd = length(this.parts);
        
        % compute total number of partitions
        dims = zeros(1, nd);
        for i = 1:nd
            dims(i) = length(this.parts{i});
        end
    end
    
    function p = getDimensions(this)
        % deprecated: use "dimensions" instead
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getDimensions'' is obsolete, use ''dimensions'' instead');
        p = partitions(this);
    end

    function dims = blockPartition(this, dim)
        % Return the dimensions of the block in the specified dimension
        %
        % DIMS = blockPartition(BD, IND)
        % DIMS is an IntegerPartition
        dims = this.parts{dim};
    end

    function dims = getBlockDimensions(this, dim)
        % Return the dimensions of the block in the specified dimension
        %
        % DIMS = getBlockDimensions(BD, IND)
        % DIMS is an IntegerPartition
        dims = this.parts{dim};
    end
    
    function dim = dimensionality(this)
        % Return the number of dimensions
        dim = length(this.parts);
    end
    
    function siz = getSize(this, varargin)
        % Return the size (number of matrix elements) in each direction
        %
        % deprecated: use size instead
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getSize'' is obsolete, use ''size'' instead');
        
        siz = size(this, varargin{:});
    end
        
    function siz = blockSize(this, varargin)
        % Return the number of size blocks of this BlockDimensions object
        %
        % BSIZ = blockSize(BD);
        % Return block size, as a row vector of block numbers in each
        % direction.
        %
        % BSIZ = blockSize(BD, DIM);
        % Return the number of blocks in the specifed direction(s).
        
        nd = length(this.parts);
        
        if isempty(varargin)
            % block size in each direction (return a vector)
            siz = zeros(1, nd);
            for i = 1:nd
                siz(i) = length(this.parts{i});
            end
            
        else
            % block size in the specified direction(s)
            % (return either a scalar, or a vector if several directions
            % are specified)
            dim = varargin{1};
            siz = zeros(1, length(dim));
            for i = 1:length(dim)
                if dim(i) > nd
                    error(sprintf(...
                        'dimension %d is too high, should be less than', dim(i), nd)); %#ok<SPERR>
                end
                siz(i) = length(this.parts{dim(i)});
            end
        end
    end
    
    function n = blockNumber(this)
        % Return the total number of blocks defined by this BlockDimensions
        
        % iterate over dimensions, to multiplies block numbers
        n = 1;
        for i = 1:length(this.parts)
            n = n * length(this.parts{i});
        end
    end
    
    function n = getBlockNumber(this, varargin)
        % Return the total number of blocks
        %
        % deprecated: use blockSize instead
        
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getBlockNumber'' is obsolete, use ''blockSize'' instead');
       
        if isempty(varargin)
            % compute total number of blocks
            n = 1;
            for i = 1:length(this.parts)
                n = n * length(this.parts{i});
            end
        else
            % returns the number of blocks only in the specified
            % dimension(s)
            dim = varargin{1};
            n = zeros(1, length(dim));
            for i = 1:length(dim)
                n(i) = length(this.parts{dim(i)});
            end
        end
    end
    
    function n = getBlockNumbers(this)
        % Return the number of blocks in each dimension
        %
        % N = getBlockNumbers(BD);
        % N is a 1-by-ND row vector
        %
        % deprecated: use blockSize instead
        
        warning('BlockMatrixToolbox:deprecated', ...
            'method ''getBlockNumbers'' is obsolete, use ''blockSize'' instead');
        
        nd = length(this.parts);
        n = zeros(1, nd);
        for i = 1:nd
            n(i) = length(this.parts{i});
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
        unif1 = isUniform(this.parts{1});
        unif2 = isUniform(this.parts{2});
        tf = unif1 && unif2;
    end
    
    function tf = isVectorBlock(this)
        % Check if a BlockMatrix is divided in 1 block in at least one direction
        tf = any(blockSize(this) == 1);
    end
end


%% overload some native methods

methods
    function siz = size(this, varargin)
        % Return the size (number of matrix elements) in each direction
        %
        % SIZ = size(BD)
        % Returns the size as a 1-by-ND row vector, where ND is the
        % dimensionality of this BlockDimensions.
        %
        % SIZ = size(BD, DIM)
        % Returns the size in the specified dimension. DIM should be an
        % integer between 1 and ND
        %

        % number of dimensions
        nd = length(this.parts);
        
        if isempty(varargin)
            % return dimension vector
            siz = zeros(1, nd);
            for i = 1:nd
                siz(i) = sum(this.parts{i});
            end
            
        else
            % return size in the specified direction
            dim = varargin{1};
            siz = zeros(1, length(dim));
            for i = 1:length(dim)
                if dim(i) > nd
                    error(sprintf(...
                        'dimension %d is too high, should be less than', dim(i), nd)); %#ok<SPERR>
                end
                siz(i) = sum(this.parts{dim(i)});
            end
        end
    end
    
    function res = transpose(this)
        % transpose the block dimensions
        res = ctranspose(this);
    end
    
    function res = ctranspose(this)
        % overload the transpose operator for BlockDimensions objects
        nd = length(this.parts);
        if nd ~= 2
            error('transpose is defined only for 2D BlockDimensions objects');
        end
        
        parts2 = cell(1, nd);
        parts2{1} = this.parts{2};
        parts2{2} = this.parts{1};
        res = BlockDimensions(parts2);
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
        
        nd = length(this.parts);
        partsH = this.parts{2};
        
        for i = 1:length(varargin)
            var = varargin{i};
            % additional BlockDimensions should have same dimensionality
            % and same block size in other dimensions
            if length(var.parts) ~= nd
                error('Other BlockDimensions should have same dimensionality');
            end
            for iDim = [1 3:nd]
                if length(this.parts{iDim}) ~= length(var.parts{iDim})
                    error('BlockDimensions should have same block number in dimension %d', iDim); 
                end
                if this.parts{iDim} ~= var.parts{iDim}
                    error('BlockDimensions should have same block sizes in dimension %d', iDim); 
                end
            end
            
            partsH = [partsH var.parts{2}]; %#ok<AGROW>
        end
        
        newParts = [this.parts(1) {partsH} this.parts(3:end)];
        res = BlockDimensions(newParts);
    end

    function res = vertcat(this, varargin)
        % Overload the vertical concatenation operator
        
        nd = length(this.parts);
        partsV = this.parts{1};
        
        for i = 1:length(varargin)
            var = varargin{i};
            % additional BlockDimensions should have same dimensionality
            % and same block size in other dimensions
            if length(var.parts) ~= nd
                error('Other BlockDimensions should have same dimensionality');
            end
            for iDim = 2:nd
                if length(this.parts{iDim}) ~= length(var.parts{iDim})
                    error('BlockDimensions should have same block number in dimension %d', iDim); 
                end
                if this.parts{iDim} ~= var.parts{iDim}
                    error('BlockDimensions should have same block sizes in dimension %d', iDim); 
                end
            end
            
            partsV = [partsV var.parts{1}]; %#ok<AGROW>
        end
        
        newParts = [{partsV} this.parts(2:end)];
        res = BlockDimensions(newParts);
    end
    
    function b = eq(this, that)
        % Compare integer partitions of two block dimensions
        %
        % BD1 = BlockDimensions( {[2 2], [2, 3, 2]} );
        % BD2 = BlockDimensions( {[2 2], [2, 2, 3]} );
        % BD1 == BD2
        % ans = 
        %      1    0
        %
        
        % both objects should be of the right class
        if ~isa(this, 'BlockDimensions') || ~isa(that, 'BlockDimensions')
            b = false;
            return;
        end
        
        % both objects should have same dimensionality
        if length(this.parts) ~= length(that.parts)
            error('BlockDimensions objects must have same dimensionality');
        end
        
        % all partitions should be equal
        b = cellfun(@eq, this.parts, that.parts);
    end
    
    function b = ne(this, that)
        % Test whether two block dimensions are different or not
        b = ~eq(this, that);
    end    
    
    function varargout = subsref(this, subs)
        % Return the integer partition for a given dimension
        %
        % DIMS = BlockDimensions({[1 2 1], [3 4]});
        % DIMS{1}
        % ans =
        % IntegerPartition object with 3 terms
        %     (1, 2, 1)
        
        
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
            
            % stop here
            return;
            
        elseif strcmp(type, '()')
            % Process parens indexing
            error('parens indexing of BlockDimensions is not supported');
            
        elseif strcmp(type, '{}')
            % Process braces indexing
            ns = length(s1.subs);
            if ns == 1
                % returns integer partition of corresponding dimension
                intPart = this.parts{s1.subs{1}};
                if length(subs) == 1
                    varargout{1} = intPart;
                else
                    % process other calls to subsref
                    varargout = cell(1);
                    [varargout{:}] = subsref(intPart, subs(2:end));
                end
                
            else
                error('Only linear indexing is allowed for BlockDimensions');
            end
        end
        
    end
    
    function n = numArgumentsFromSubscript(this,~,~)
        % Need to overload this to allow proper braces indexing
        n = numel(this);
    end
end

%% Display methods

methods
    function disp(this)
        % Display the content of this BlockDimensions object
        
        % loose format: display more empty lines
        isLoose = strcmp(get(0, 'FormatSpacing'), 'loose');
        
        % get dimensionality
        nd = dimensionality(this);
        
        % Display information on block sizes in each dimension
        disp(sprintf('BlockDimensions object with %d dimensions', nd)); %#ok<DSPS>
        disp(char(this));
        
        if isLoose
            fprintf('\n');
        end
    end

    function buffer = char(this)
        % Convert BlockDimension object to string representation
        
        nd = length(this.parts);
        buffer = ['( ' char(this.parts{1})];
        for i = 2:nd
            buffer = [buffer ', ' char(this.parts{i})]; %#ok<AGROW>
        end
        buffer = [buffer ' )'];
    end
end % end methods

end % end classdef

