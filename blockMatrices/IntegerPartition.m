classdef IntegerPartition < handle
%INTEGERPARTITION Store an ordered partition of an integer
%
%   Class IntegerPartition
%   Store an ordered partition of an integer, as a list of integer terms.
%   Terms should be positive integers.
%
%   Example
%     IP = IntegerPartition([2, 3, 4]);
%     length(IP)
%     ans =
%         8
%     term(IP, 3)
%     ans = 
%         4
%
%   See also
%     BlockDimensions

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-04-10,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % the partition of the integer
    % given as a 1-by-N row vector of partitions
    terms;
    
end % end properties


%% Static methods
methods (Static)
    function res = ones(n)
        % Create a new partition formed only by ones
        %
        % Example
        % I4 = IntegerPartition.ones(4)
        % I4 =
        % IntegerPartition object with 4 terms
        %     (1, 1, 1, 1)
        res = IntegerPartition(ones(1, n));
    end
    
    function count = countPartitions(n, s)
        % Count the number of partitions of integer N into S parts.
        %
        % NP = IntegerPartition.countPartitions(N, S)
        % N:    the integer to be partitioned
        % S:    the number of partitions
        % NP:   the total number of ways the integer N can be partitioned
        %       into S parts.
        %
        % Example
        %   IntegerPartition.countPartitions(10, 4)
        %   ans =
        %       84
        %
        % See also
        %   randomPartition, choosePartition
        
        % get the singleton instance of the dictionary
        dict = IntegerPartition.getDictionaryInstance();
        
        % compute key for current case
        key = sprintf('%d.%d', n, s);
        
        % check if result was already computed
        if isKey(dict, key)
            count = dict(key);
            return;
        end
        
        % case of S partitions of 1 element, or 1 partition -> only one possibility
        if n == s || s == 1
            count = 1;
            dict(key) = count; %#ok<NASGU>
            return;
        end
        
        count = 0;
        for k = 1:(n-s+1)
            count = count + IntegerPartition.countPartitions(n-k, s-1);
        end
        dict(key) = count; %#ok<NASGU>
    end
    
    function ip = choosePartition(n, s, r)
        %CHOOSEPARTITION Choose a partition from its index
        %
        %   C = IntegerPartition.choosePartition(N, S, R)
        %   N: the integer to be partitioned
        %   S: the number of partitions
        %   R: the index of the partition, between 1 and countPartitions(N, S)
        %
        %   Example
        %     IntegerPartition.choosePartition(10, 4, 12)
        %
        %   See also
        %     countPartitions, randomPartition
        
        % get the singleton instance of the dictionary
        dict = IntegerPartition.getDictionaryInstance();
        
        % allocate memory for result
        parts = zeros(1, s);
        
        % init for first iteration
        ni = n;
        si = s;
        ri = r;
        
        % iterate over the partitions
        for i = 1:s-1
            % compute maximal admissible value for left-most partition
            kiMax = ni + 1 - si;
            
            % for each admissible value, retrieve the number of partitions
            partsi = zeros(1, kiMax);
            for ki = 1:kiMax
                % compute key for current case
                key = sprintf('%d.%d', ni-ki, si-1);
                if isKey(dict, key)
                    partsi(ki) = dict(key);
                else
                    disp(sprintf('strange, key %s is not stored...', key)); %#ok<DSPS>
                    partsi(ki) = IntegerPartition.countPartitions(ni, si);
                    dict(key) = partsi(ki);
                end
            end
            
            % identify the value of current partition corresponding to chosen
            % integer
            ind = find(cumsum(partsi) >= ri, 1, 'first');
            parts(i) = ind;
            
            % prepare for next iteration
            ni = ni - ind;
            si = si - 1;
            ri = ri - sum(partsi(1:ind-1));
        end
        
        % compute the value of the last partition
        parts(end) = n - sum(parts(1:end-1));
        
        % convert to IntegerPartition class
        ip = IntegerPartition(parts);
    end
    
    function ip = randomPartition(n, s)
        %RANDOMPARTITION Choose a random partition of integer N into S parts
        %
        %   PARTS = IntegerPartition.randomPartition(N, S)
        %   N: the integer to be partitioned
        %   S: the number of partitions
        %
        %   Example
        %     IntegerPartition.randomPartition(10, 4)
        %
        %   See also
        %     countPartitions, choosePartition
        
        np = IntegerPartition.countPartitions(n, s);
        index = randi(np);
        ip = IntegerPartition.choosePartition(n, s, index);
    end
end

% methods (Static, Access = private)
methods (Static)
    function dict = getDictionaryInstance()
        % Return the instance of dictionary storing partition counts
        dict = getDictionary(PartitionCountDictionary.getInstance());
    end
end

%% Constructor
methods
    function this = IntegerPartition(varargin)
    % Constructor for IntegerPartition class
    % 
    % IP = IntegerPartition(TERMS)
    % where TERMS is a row vector of positive integers, initialize the
    % partition with the given terms.
    %
    % IP = IntegerPartition(IP0)
    % Copy constructor
    
        if nargin == 1 && isnumeric(varargin{1})
            % initialisation constructor from numeric array
            var1 = varargin{1};
            if any(var1 <= 0) || any(round(var1) ~= var1)
                error('Requires only positive integers');
            end
            
            % ensures terms are stored in row vector
            this.terms = varargin{1}(:)';
            
        elseif nargin == 1 && isa(varargin{1}, 'IntegerPartition')
            % copy constructor
            var1 = varargin{1};
            this.terms = var1.terms;
            
        else
            error('Requires an initialisation array');
        end

    end

end % end constructors


%% Methods
methods
    function p = term(this, index)
        % Return the size of the i-th term of the partition
        %
        % Example
        %   IP = IntegerPartition([3 4 5])
        %   term(IP, 3)
        %   ans = 
        %       3
        %
        % see also
        %   length
        
        p = this.terms(index);
    end
    
    function n = integer(this)
        % Return the value of the partitioned integer
        % deprecated: use sum instead
        warning('deprecated: use sum instead');
        n = sum(this.terms);
    end
    
    function inds = blockIndices(this, blockInds)
        % Return the linear indices of the elements in the i-th block
        %
        % EINDS = blockIndices(IP, BINDS)
        % IP is an integer partition, BINDS is the array of block indices,
        % and EINDS is the array of element indices referenced by BINDS.
        %
        %
        % Example
        %   IP = IntegerPartition([2 3 4 2])
        %   blockIndices(IP, 3)
        %   ans = 
        %       6   7   8   9
        %
        %   blockIndices(IP, [1 4])
        %   ans =
        %       1   2  10  11
        %
        % see also
        %   term, length
        
        inds = [];
        for i = 1:length(blockInds)
            index = blockInds(i);
            inds_i = (1:this.terms(index)) + sum(this.terms(1:index-1));
            inds = [inds inds_i]; %#ok<AGROW>
        end
    end

end % end methods

%% boolean methods to identify the type of partition
methods
    function tf = isUniform(this)
        % Return true if all terms are equal
        tf = all(this.terms == this.terms(1));
    end
    
    function tf = isScalar(this)
        % Return true if the length of the partition equals one
        tf = length(this.terms) == 1;
    end
    
    function tf = isOnes(this)
        % Return true if all terms equal one
        % (the method isUniform will return true as well). 
        tf = all(this.terms == 1);
    end
end % end methods

%% Overload some native methods

methods
    function n = length(this)
        % Return the number of terms of this partition
        n = length(this.terms);
    end
    
    function n = sum(this)
        % Return the sum of the terms
        % (returns the same result as the "integer" function)
        n = sum(this.terms);
    end

    function res = mtimes(this, that)
        % Multiply a partition by a scalar integer
        %
        % P2 = mtimes(P1, S)
        % P2 = mtimes(S, P1)
        %
        % Example
        % P = IntegerPartition([1 2 3]);
        % P2 = P * 3
        % P2 = 
        % IntegerPartition object with 3 terms
        %   ( 3, 6, 9)
        
        
        % one of the two arguments is an integer
        % -> identify arguments
        if isa(this, 'IntegerPartition')
            part = this;
            arg = that;
        else
            part = that;
            arg = this;
        end
        
        % check validity (scalar and integer)
        if ~isscalar(arg) || mod(arg, 1) ~= 0
            error('second argument must be a scalar integer');
        end
        
        % create result
        newTerms = part.terms * arg;
        res = IntegerPartition(newTerms);
    end
    
    function res = times(this, that)
        % Multiply two partitions element-wise
        %
        % P3 = times(P1, P2)
        % P3 = P1 .* P2
        
        if ~isa(this, 'IntegerPartition') || ~isa(that, 'IntegerPartition')
            error('Both arguments must be IntegerPartition');
        end
        
        % Both arguments are instances of integer partition
        % -> use element-wise multiplication
        
        % check length
        if length(this) ~= length(that)
            error('The two partitions must have the same length');
        end
        
        newTerms = this.terms .* that.terms;
        res = IntegerPartition(newTerms);
    end
    
    function res = plus(this, that)
        % Add two partitions element-wise
        %
        % P3 = plus(P1, P2)
        % P3 = P1 + P2
        
        if isa(this, 'IntegerPartition') && isa(that, 'IntegerPartition')
            % Both arguments are instances of integer partition
            % -> use element-wise addition

            % check length
            if length(this) ~= length(that)
                error('The two partitions must have the same length');
            end

            newTerms = this.terms + that.terms;
            res = IntegerPartition(newTerms);
        else
            
            % one of the two arguments is numeric
            % -> identify arguments
            if isa(this, 'IntegerPartition')
                part = this;
                arg = that;
            else
                part = that;
                arg = this;
            end
            
            newTerms = part.terms + arg;
            res = IntegerPartition(newTerms);
        end
    end
    
    function res = mrdivide(this, arg)
        % Divide partiton terms by an integer
        
        if ~isscalar(arg) || mod(arg, 1) ~= 0
            error('second argument must be an integer');
        end
        if any(mod(this.terms, arg) ~= 0)
            error('at least one term is not divisible by %d', arg);
        end
        
        newTerms = this.terms / arg;
        res = IntegerPartition(newTerms);
    end
    
    function res = horzcat(this, varargin)
        % Overload the horizontal concatenation operator
        
        newTerms = this.terms;
        
        for i = 1:length(varargin)
            that = varargin{i};
            if ~isa(that, 'IntegerPartition')
                error(['Additional argument should be an IntegerPartition, not a ' classname(that)]);
            end
            newTerms = [newTerms that.terms]; %#ok<AGROW>
        end
        
        res = IntegerPartition(newTerms);
    end
    
    function varargout = subsref(this, subs)
        % Return the term of this partition at the given index
        %
        % P = IntegerPartition([2 3 2]);
        % P(2)
        % ans =
        %     3
        
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
            
            varargout{1} = 0;
            
            % check number of indices
            ns = length(s1.subs);
            if ns == 1
                % returns the requested terms as a row vector
                varargout{1} = this.terms(s1.subs{1});
            else
                error('Only linear indexing is allowed for IntegerPartition');
            end
            
        else
            % Process braces indexing
            
            ns = length(s1.subs);
            if ns == 1
                % returns the requested terms as a new IntegerPartition
                newTerms = this.terms(s1.subs{1});
                varargout = {IntegerPartition(newTerms)};
            else
                error('Only linear indexing is allowed for IntegerPartition');
            end
        end
            
    end
    
    function n = numArgumentsFromSubscript(this,~,~)
        % Need to overload this to allow proper braces indexing
        n = numel(this);
    end

    function b = eq(this, that)
        % Test whether two compositions are the same or not
        
        if ~isa(this, 'IntegerPartition') || ~isa(that, 'IntegerPartition')
            b = false;
            return;
        end
        
        if length(this.terms) ~= length(that.terms)
            b = false;
            return;
        end
        
        b = all(this.terms == that.terms);
    end
    
    function b = ne(this, that)
        % tests whether two compositions are different or not
        b = ~eq(this, that);
    end
    
end


%% Display methods

methods
    function disp(this)
        nd = length(this.terms);
        disp(sprintf('IntegerPartition object with %d terms', nd)); %#ok<DSPS>
        disp(['    ' char(this)]);
    end
    
    function buffer = char(this)
        % convert to string representation
        
        n = length(this.terms);
        pattern = ['(%d' repmat(', %d', 1, n-1) ')'];
        buffer = sprintf(pattern, this.terms);
    end

end

end % end classdef

