classdef Iteration < handle
%ITERATION Performs iteration of block-matrix factorization algorithm
%
%   Class Iteration
%
%   Example
%   Iteration
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-11-27,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % the working block-matrix
    data;
    
    % the objective function to minimize, as a function handle
    criterion;
    
    % the initialisation vector, as an instance of BlockMatrix
    init;
    
    % the resulting vector, as an instance of BlockMatrix
    solution;
    
    
end % end properties


%% Constructor
methods
    function this = Iteration(varargin)
    % Constructor for Iteration class

    end

end % end constructors


%% Methods
methods
end % end methods

end % end classdef

