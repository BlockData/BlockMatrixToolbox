classdef JacobiBlockPower < BlockPowerAlgo
%JACOBIBLOCKPOWER Jacobi algorithm for solving block power algorithms
%
%   STATE = JacobiBlockPower(BM)
%   Creates a new instance of Jacobi Block Power iteration algorithm, using
%   the specified Block-Matrix BM for representing the problem, and an
%   optional Block-Vector representing the initial state of the algorithm.
%
%   The algorithm can be used that way:
%   NEWSTATE = iterate(STATE)
%
%   Example
%     data0 = round(rand(4, 7) * 1000) / 1000;
%     data = BlockMatrix.create(data0, 4, [2 3 2]);
%     AA = blockProduct_uu(data', data);
%     AA.data = max(AA.data, AA.data'); % ensure symmetry of the matrix
%     algo = JacobiBlockPower(AA);
%
%
%   See also
%     GaussBlockPower

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-01-15,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - BIA-BIBS.


%% Properties
properties
    % the BlockMatrix representing the problem
    data;

    % a function handle for computing matrix from A and u. Default is A.
    core;

    % the type of block product to apply on (block)matrix and (block)vector
    % default is "uu"
    productType = 'uu';

    % a function handle that computes new value of vector from the product
    % A *_(w1,w2) u.
    updateFunction;

    % the type of norm used to normalize vector. Default is L2 norm.
    normType = 2;

end % end properties


%% Constructor
methods
    function this = JacobiBlockPower(A, varargin)
        % Constructor for JacobiBlockPower class
        %
        % Usage:
        %   STATE = JacobiBlockPower(MAT);
        %   MAT is a N-by-P BlockMatrix.
        %

        % Empty constructor, to allow creation of arrays
        if nargin == 0
            return;
        end

        % Copy constructor
        if isa(A, 'JacobiBlockPower')
            this.data   = A.data;
            this.core   = A.core;

            this.productType    = A.productType;
            this.updateFunction = A.updateFunction;
            this.normType       = A.normType;

            return;
        end

        % check type of first argument
        if ~isa(A, 'BlockMatrix')
            error('First argument should be a block-matrix');
        end

        % check matrix validity.
        % PositiveDefinite Matrices are not tested for now.
        [n1, n2] = size(A);
        if n1 ~= n2
            error('Requires a square matrix');
        end
        if ~isSymmetric(A)
            warning('Requires a symmetric matrix');
        end
        if ~isPositiveDefinite(A)
            warning('Requires a positive definite matrix');
        end
        this.data = A;

        % default core function simply returns the original matrix.
        this.core = @(A,u) A;

        % create the default update function
        this.updateFunction = @(Au) blockProduct_hs( 1 ./ blockNorm(Au), Au );
    end

end % end constructors


%% Algorithm monitoring methods
methods

    function stateList = convergence(this, state, varargin)
        % Iterates this algorithm until a stopping criterion is found
        %
        % PATH = convergence(ALGO, U0);
        %
        % PATH = convergence(..., PNAME, PVALUE)
        % Specified one or several optional parameter as name-value pairs.
        %
        % List of available parameters
        % * maxIterNumber:  the maximum number of iteration (default 100)
        % * residTol:       the tolerance on residuals, as the norm of the
        %       block-norm of the difference between two successive vectors.
        %       Default value is 1e-8.
        % * criteriumTol:   the tolerance on the difference between two
        %       successive values of the computed eigen value. Default
        %       value is 1 e-8.
        %
        % See also
        %   solve, blockPowerOptions

        % parse optimization options
        options = blockPowerOptions(varargin{:});

        % initialize list of state
        stateList = {state};

        % iterate until residual is acceptable
        for iIter = 1:options.maxIterNumber
            % performs one iteration, and agglomerate
            state = this.iterate(state);
            stateList = [stateList {state}]; %#ok<AGROW>

            % test the tolerance on residual
            if state.residual < options.residTol
                fprintf('converged to residual after %d iteration(s)\n', iIter);
                return;
            end

            %             % test the tolerance on eigen value
            %             if eigenValue(state) < options.eigenTol
            %                 fprintf('converged to eigenValue after %d iteration(s)\n', iIter);
            %                 return;
            %             end
        end

        fprintf('Reached maximum number of iterations (%d)\n', iIter);
    end

    function t = tolerance(this, state) %#ok<INUSL>
        % Returns the residual, or difference with previous iteration
        t = state.residual;
    end

    function varargout = monotony(this, varargin)
        % Computes monotony of current state using: u' * A(u) * u
        %
        % Usage:
        %    monotony(ALGO)
        %    monotony(..., 'nIter', N);
        %    Specifies the number of iterations
        %    monotony(..., 'display', DISP)
        %    Specifies whether the resulting curve should be displayed
        %    (DISP = 'on') or not (DISP = 'off').
        %

        % default value for iteration number
        nIter = 100;
        display = 'on';

        % get iteration number if given as numeric
        if ~isempty(varargin) && isnumeric(varargin{1})
            nIter = varargin{1};
            varargin(1) = [];
        end

        % parses input argument
        while length(varargin) > 1
            paramName = varargin{1};
            if strcmpi(paramName, 'nIter')
                nIter = varargin{2};
            elseif strcmpi(paramName, 'display')
                display = varargin{2};
            else
                error(['Unknown parameter name: ' paramName]);
            end
            varargin(1:2) = [];
        end

        % initialize reuslt array
        values = zeros(1, nIter);

        % iterate algorithm for the given number of iterations
        state = this;
        for i = 1:nIter
            u = state.vector;
            Au = computeProduct(state, u);
            values(i) = getMatrix(u' * Au);
            state = state.next();
        end

        % eventually displays the result
        if strcmp(display, 'on')
            figure;
            plot(values);
            xlabel('Iterations');
        end

        % eventually returns the computed values
        if nargout > 0
            varargout{1} = values;
        end
    end

    function bestSolution = globality(this, varargin)
        % Performs several trials to fin the best solution in unit hypercube
        %
        % usage
        %    STATE = globality(ALGO);
        %
        %    STATE = globality(ALGO, NSIMS);
        %    Where NSIMS is the number of simulations.
        %
        %    STATE = globality(..., OPT_NAME, OPT_VAL);
        %    Specifies optional parameters as pairs of name-value. See
        %    blockPowerOptions for details.
        %
        % see also
        %    isGlobal

        % parse number of simulations
        nSimuls = 10;
        if ~isempty(varargin) && isnumeric(varargin{1})
            nSimuls = varargin{1};
            varargin(1) = [];
        end

        % integer partition of the vector
        vparts = blockDimensions(this.data, 2);

        % block dimensions of the vector
        vdims = BlockDimensions({vparts, IntegerPartition(1)});

        % init residual to max value to force keeping first solution
        bestResidual = Inf;

        % run the simulations
        for i = 1:nSimuls
            % compute a new solution
            vector = BlockMatrix.create(rand(size(vdims)), vdims);
            state = BlockPowerAlgoState(vector);
            solution = solve(this, state, varargin{:});

            % compare residuals, and keep best solution
            if solution.residual < bestResidual
                bestSolution = solution;
                bestResidual = solution.residual;
            end
        end
    end
end


%% Iteration methods
methods
    function state = solve(this, state, varargin)
        % Iterates this algorithm until a stopping criterion is found
        %
        % U = solve(ALGO, U0);
        % ALGO is the BlockPowerAlgo instance
        % U0 is an instance of BlockPowerAlgoState, or a BlockVector.
        %
        % U = solve(..., PNAME, PVALUE)
        % Specified one or several optional parameter as name-value pairs.
        %
        % List of available parameters
        % * maxIterNumber:  the maximum number of iteration (default 100)
        % * residTol:       the tolerance on residuals, as the norm of the
        %       block-norm of the difference between two successive vectors.
        %       Default value is 1e-8.
        % * eigenTol:       the tolerance on the difference between two
        %       successive values of the computed eigen value. Default
        %       value is 1 e-8.
        %
        % See also
        %   convergence, blockPowerOptions

        % ensure first argument is AlgoState
        if ~isa(state, 'BlockPowerAlgoState')
            % try to convert argument to AlgoState instance
            state = BlockPowerAlgoState(state);
        end

        % parse optimization options
        options = blockPowerOptions(varargin{:});

        % iterate until residual is acceptable
        for iIter = 1:options.maxIterNumber
            % performs one iteration, and get residual
            state = this.iterate(state);

            % test the tolerance on residual
            if state.residual < options.residTol
                fprintf('converged to residual after %d iteration(s)\n', iIter);
                return;
            end

            %             % test the tolerance on eigen value
            %             if eigenValue(state) < options.eigenTol
            %                 fprintf('converged to eigenValue after %d iteration(s)\n', iIter);
            %                 return;
            %             end
        end

        fprintf('Reached maximum number of iterations (%d)\n', iIter);
    end

    function newState = iterate(this, state)
        % Performs a single iteration of the (Block-)Power Algorithm
        %
        % NEWSTATE = iterate(ALGO, STATE)
        % where STATE is a correctly initialized JacobiBlockPower
        % algorithm, returns the new state of the algorithm, as an instance
        % of BlockPowerAlgoState.
        %

        % extract vector
        qq = state.vector;

        % compute the matrix from the core function and the input data
        A = this.core(this.data, qq);

        % performs block-product on current vector
        q = blockProduct(A, qq, this.productType);

        % block normalization
        q = this.updateFunction(q);
        % usually:
        % q = blockProduct_hs(1./blockNorm(q), q);

        % create algorithm state data structure
        newState = BlockPowerAlgoState(q);

        % compute residual
        resid = norm(blockNorm(q - qq), this.normType);
        newState.residual = resid;
    end
end


%% Utility methods
methods
    function lambdas = pseudoEigenValues(this, state)
        % Returns the pseudo eigen values of current state
        %
        %   Usage
        %   LAMBDAS = pseudoEigenValues(ALGO, STATE)
        %
        % (renamed from "function lambdas = stationarity(this)")

        u = state.vector;
        Au = computeProduct(this, u);
        blockLambdas = blockNorm(Au) ./ blockNorm(u);
        lambdas = getMatrix(blockLambdas);
    end

    function res = isGlobal(this, state)
        % Returns 1 of eig(A - lambda * I) or NaN otherwise.
        %
        %   RES = isGlobal(ALGO, STATE)
        %
        % see also
        %    globality

        % get vector and core matrix
        u = state.vector;
        A = coreMatrix(this, u);

        % compute lambda for current iteration
        Au = computeProduct(this, u);
        blockLambdas = blockNorm(Au) ./ blockNorm(u);
        lambdas = getMatrix(blockLambdas);

        % extract block-dimension
        dims = blockDimensions(A);
        dim1 = dims{1};

        % initialize the lambda_k * I
        nBlocks = length(dim1);
        blocks = cell(1, nBlocks);
        for i = 1:nBlocks
            blocks{i} = ones(dim1(i)) * lambdas(i);
        end

        lambdaI = BlockDiagonal(blocks);
        t = eig(getMatrix(A - lambdaI));

        if t < 0
            res = 1;
        else
            res = NaN;
        end
    end

    function lambda = eigenValue(this, state)
        % Computes the current eigen value
        %
        %   LAMBDA = eigenValues(ALGO, STATE)

        q = blockProduct(this.data, this.vector, this.productType);
        lambda = norm(q, this.normType);
    end

    function Au = computeProduct(this, u)
        % Computes the (w1,w2)-product of core matrix by the specified vector
        A = this.core(this.data, u);
        Au = blockProduct(A, u, this.productType);
    end

    function A = coreMatrix(this, u)
        % Computes the current core matrix, from data matrix and vector U
        A = this.core(this.data, u);
    end

    function un = normalizeVector(this, u)
        % Computes normalized vector, using inner settings
        un = this.updateFunction(u);
    end
end % end methods

end % end classdef

