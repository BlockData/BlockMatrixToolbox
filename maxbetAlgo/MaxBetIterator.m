classdef MaxBetIterator < handle
%MAXBETITERATOR Iterate the MaxBet algorithm
%
%   Class MaxBetIterator.
%   This class is used by the 'maxbet_procedure3' method.
%
%   Example
%   MaxBetIterator
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
    % the data of the problem.
    % For MAXBET algorithm, corresponds to: blockProduct_uu(data', data);
    data;
    
    % the current vector. 
    % Iterations should converge towards stable value.
    vector;
    
end % end properties


%% Constructor
methods
    function this = MaxBetIterator(data, init)
    % Constructor for MaxBetIterator class
        this.data = data;
        this.vector = init;
    end

end % end constructors


%% Methods
methods
    function [q, resid] = iterate(this)
        % Perform a single MaxBet iteration
        %
        % [Q, RESID] = iterate(MBI)
        % where MBI is a correctly initialized MaxBet iterator, returns the
        % new value Q of the vector (as a BlockMatrix), and the value of
        % the residual (as a scalar).
        %
        
        % performs block-product on current vector (uses 'uu' product)
        qq = this.vector;
        q = blockProduct_uu(this.data, qq);
        
        % block normalization
        q = blockProduct_hs(1./blockNorm(q), q); 
        
        % compute residual
        resid = norm(blockNorm(q) - blockNorm(qq)); 
        
        % keep result for next iteration
        this.vector = q;
    end
    
end % end methods

end % end classdef

