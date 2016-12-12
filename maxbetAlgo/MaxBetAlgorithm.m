classdef MaxBetAlgorithm < handle
%MAXBETALGORITHM Maxbet factorization algorithm
%
%   Class MaxBetAlgorithm
%
%   Example
%   ALGO = MaxBetAlgorithm('tolerance', .01, 'maxIterNumber'
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-12-14,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    tolerance = .01;
    
    nbMaxIter = 100;
    
end % end properties


%% Constructor
methods
    function this = MaxBetAlgorithm(varargin)
        % Constructor for MaxBetAlgorithm class
        
        while length(varargin) > 1
            paramName = varargin{1};
            switch paramName
                case 'tolerance'
                    this.tolerance = varargin{2};
                case 'nbmaxiter'
                    this.nbMaxIter = varargin{2};
                otherwise
                    error(['Unknown parameter name: ' paramName]);
            end
            
            varargin(1:2) = [];
        end
    end

end % end constructors


%% Methods
methods
    function [q, resid, nIter] = factorize(this, data, q0)
        % Apply MaxBet Algorihtm to input data
        
        % compute the block-matrix corresponding to maxbet algorithm
        AA = blockProduct_uu(data', data);
        
        % create new BlockMatrix representing the normalized input vectors
        qq = blockProduct_hs(1./blockNorm(q0), q0);
        
        % create the algorithm iterator
        iterator = MaxBetIterator(AA, qq);
        
        % init residual
        resid = 1;
        
        % iterate until residual is acceptable
        nIter = 0;
        while true
            % increment iteration
            nIter = nIter + 1;
            
            % performs one iteration, and get residual
            [q, resid] = iterator.iterate();
            
            % check stopping conditions
            if resid < this.tolerance
                break;
            end
            if nIter >= this.nbMaxIter
                break;
            end
        end

    end
end % end methods

end % end classdef

