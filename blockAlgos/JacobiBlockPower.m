classdef JacobiBlockPower < BlockPowerAlgoState
%JACOBIBLOCKPOWER Jacobi algorithm for solving block power algorithms
%
%   STATE = JacobiBlockPower(BM)
%   STATE = JacobiBlockPower(BM, U0)
%   Creates a new instance of Jacobi Block Power iteration algorithm, using
%   the specified Block-Matrix BM for representing the problem, and an
%   optional Block-Vector representing the initial state of the algorithm.
%   If U0 is not specified, a block vector containing only 1 is used.
%
%   The algorithm can be used that way:
%   NEWSTATE = next(STATE)
%
%   Example
%
%
%   See also
%     GaussBlockPower

% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2016-01-15,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - BIA-BIBS.


%% Properties
properties
    % the BlockMatrix representing the problem
    data;
    
    % the block vector representing the initial solution
    vector;

    % a function handle for computing matrix from A and u. Default is A.
    core;
    
    % the residual from previous state. NaN for first iteration.
    residual = NaN;

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
        %   STATE = JacobiBlockPower(MAT, U0);
        %   MAT is a N-by-P BlockMatrix. U0 is an optional N-by-1
        %   BlockVector.
        %
        
        % Empty constructor, to allow creation of arrays
        if nargin == 0
            return;
        end
        
        % Copy constructor
        if isa(A, 'JacobiBlockPower')
            this.data   = A.data;
            this.vector = A.vector;
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
        if n1 ~= n2 || ~isSymmetric(A) % || ~isPositiveDefinite(A)
            error('Requires a symmetric positive definite matrix');
        end
        if ~isPositiveDefinite(A)
            warning('Requires a symmetric positive definite matrix');
        end
        this.data = A;
        
        % Check if initialisation vector is precised
        if nargin > 1
            % use second input argument for initial vector
            var1 = varargin{1};
            if ~isa(var1, 'BlockMatrix')
                error('Second argument should be a block-matrix');
            end
            if blockSize(var1, 2) ~= 1
                error('Second argument should be a block-matrix with one block-column');
            end
            this.vector = varargin{1};
            
        else
            % create initial vector from matrix size
            n = size(A, 1);
            vdim = blockDimensions(A, 2);
            this.vector = BlockMatrix(ones(n, 1), vdim, 1);
        end
        
        % default core function simply returns the original matrix.
        this.core = @(A,u) A;
        
        % create the default update function
        this.updateFunction = @(Au) blockProduct_hs( 1 ./ blockNorm(Au), Au );
    end

end % end constructors


%% Algorithm monitoring methods
methods
    function state = solve(this, varargin)
        % Iterates this algorithm until a stopping criterion is found
        %
        % U = solve(ALGO, U0);
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
        
        % uses first argument if this is a block-vector
        if nargin > 1 && isa(varargin{1}, 'AbstractBlockMatrix')
            this.vector = varargin{1};
            varargin(1) = [];
        end
        
        % parse optimization options
        options = blockPowerOptions(varargin{:});
        
        % initialize the state
        state = this;
        
        % iterate until residual is acceptable
        for iIter = 1:options.maxIterNumber
            % performs one iteration, and get residual
            state = state.next();
            
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

    function stateList = convergence(this, varargin)
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
        % * eigenTol:       the tolerance on the difference between two
        %       successive values of the computed eigen value. Default
        %       value is 1 e-8. 
        %
        % See also
        %   solve, blockPowerOptions
        
        % uses first argument if this is a block-vector
        if nargin > 1 && isa(varargin{1}, 'AbstractBlockMatrix')
            this.vector = varargin{1};
            varargin(1) = [];
        end
        
        % parse optimization options
        options = blockPowerOptions(varargin{:});
        
        % initialize list of state
        state = this;
        stateList = {state};
        
        % iterate until residual is acceptable
        for iIter = 1:options.maxIterNumber
            % performs one iteration, and agglomerate
            state = state.next();
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
    
    function t = tolerance(this)
        % return the residual, or difference with previous iteration
        t = this.residual;
    end
    
   function t = monotony(this)
       % compute monotony of current state using: u' * A(u) * u
       % TODO: NITER
       u = this.vector;
       Au = computeProduct(this, u);
       t = u' * Au;
   end
   
   function lambdas = stationarity(this)
       % return the eigen values of current state
       u = this.vector;
       Au = computeProduct(this, u);
       blockLambdas = blockNorm(Au) ./ blockNorm(u);
       lambdas = getMatrix(blockLambdas);       
   end
   
   function t = globality(this)
       % return eig(A - lambda * I)
       
       % get vector and core matrix
       u = this.vector;
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
   end
end


%% Iteration methods
methods  
    function res = next(this)
        % Performs a single iteration of the (Block-)Power Algorithm
        %
        % NEWSTATE = iterate(STATE)
        % where STATE is a correctly initialized JacobiBlockPower
        % algorithm, returns the new state of the algorithm, as an instance
        % of JacobiPowerBlock. 
        %
        
        % extract vector
        qq = this.vector;

        % compute the matrix from the core function and the input data
        A = this.core(this.data, qq);
        
        % performs block-product on current vector
        q = blockProduct(A, qq, this.productType);
        
        % block normalization
        q = this.updateFunction(q);
        % usually:
        % q = blockProduct_hs(1./blockNorm(q), q); 
        
        % create algorithm state data structure
        res = JacobiBlockPower(this);
        res.vector = q;

        % compute residual
        resid = norm(blockNorm(q - qq), this.normType); 
        res.residual = resid;
        
    end
end

%% Utility methods
methods
    function lambda = eigenValue(this)
        % Compute the current eigen value
        q = blockProduct(this.data, this.vector, this.productType);
        lambda = norm(q, this.normType);
    end
    
    function Au = computeProduct(this, u)
        % Compute the (w1,w2)-product of core matrix by the specified vector
        A = this.core(this.data, u);
        Au = blockProduct(A, u, this.productType);
    end
    
    function A = coreMatrix(this, u)
        % Compute the current core matrix, from data matrix and vector U 
        A = this.core(this.data, u);
    end
    
    function un = normalizeVector(this, u)
        % Compute normalized vector, using inner settings
        un = this.updateFunction(u);
    end
    
    
    %% Methods for studying the algorithm 
    
%     function monotony(this, v0, varargin)
%         % Display monotony of this algorithm
%         %
%         % usage:
%         %   monotony(ALGO, V0);
%         % where ALGO if the instance of Block power algorithm, and V0 is
%         % the initial value of the vector
%         %
%         % Example
%         %   % create block matrix for problem
%         %   mdims = BlockDimensions({40, [20 30 20]});
%         %   data = BlockMatrix(rand(40, 70), mdims);
%         %   AA = blockProduct_uu(data', data);
%         %   % create block matrix for initial vector
%         %   vdims = BlockDimensions({[20 30 20], 1});
%         %   q0 = BlockMatrix(rand(70, 1), vdims);
%         %   qq = blockProduct_hs(1./blockNorm(q0), q0);
%         %   % compute algorithm monotony
%         %   algo = JacobiBlockPower(AA, qq);
%         %   monotony(algo, qq);
%         
%         % parse optimization options
%         options = blockPowerOptions(varargin{:});
%         
%         % initialize algo
%         this.vector = v0;
%         
%         % initialize display
%         nIter = options.maxIterNumber;
%         lambdaList = zeros(nIter, 1);
%         
%         % iterate
%         for iIter = 1:nIter
%             iterate(this);
%             lambdaList(iIter) = eigenValue(this);
%         end
%         
%         
%         % display result
%         figure; set(gca, 'fontsize', 14); hold on;
%         plot([1 nIter], lambdaList([end end]), 'k');
%         plot(lambdaList, 'color', 'b', 'linewidth', 2);
%         xlabel('Iteration Number');
%         ylabel('Eigen value estimation');
% 
%     end
%     function stationarity(this, v0, varargin)
%         % Display stationarity of this algorithm
%         %
%         % usage:
%         %   stationarity(ALGO, V0);
%         % where ALGO if the instance of Block power algorithm, and V0 is
%         % the initial value of the vector
%         %
%         % Example
%         %   % create block matrix for problem
%         %   mdims = BlockDimensions({40, [20 30 20]});
%         %   data = BlockMatrix(rand(40, 70), mdims);
%         %   AA = blockProduct_uu(data', data);
%         %   % create block matrix for initial vector
%         %   vdims = BlockDimensions({[20 30 20], 1});
%         %   q0 = BlockMatrix(rand(70, 1), vdims);
%         %   qq = blockProduct_hs(1./blockNorm(q0), q0);
%         %   % compute algorithm monotony
%         %   algo = JacobiBlockPower(AA, qq);
%         %   stationarity(algo, qq);
%         
%         % parse optimization options
%         options = blockPowerOptions(varargin{:});
%         
%         % initialize algo
%         this.vector = v0;
%         
%         % initialize display
%         nIter = options.maxIterNumber;
%         normList = zeros(nIter, 1);
%         
%         % iterate
%         for iIter = 1:nIter
%             iterate(this);
%             
%             u = this.vector;
%             u2 = computeProduct(this, u);
%             normList(iIter) = norm(blockNorm(u2));
%         end
%         
%         
%         % display result
%         figure; set(gca, 'fontsize', 14); hold on;
%         plot([1 nIter], normList([end end]), 'k');
%         plot(normList, 'color', 'b', 'linewidth', 2);
%         xlabel('Iteration Number');
%         ylabel('Eigen value estimation');
% 
%     end
    
end % end methods

end % end classdef

