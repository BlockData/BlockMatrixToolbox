classdef JacobiCoreBlockPower < BlockPowerAlgo
%JACOBICOREBLOCKPOWER Jacobi algorithm with kernel for solving block power algorithms
%
%   Class JacobiCoreBlockPower
%
%   Example
%   JacobiCoreBlockPower
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-01-15,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - BIA-BIBS.


%% Properties
properties
    % the BlockCore object representing the kernel
    core;
    
    % the block vector representing the initial solution
    init;
    
    % the current solution
    solution;
    
end % end properties


%% Constructor
methods
    function this = JacobiCoreBlockPower(varargin)
    % Constructor for JacobiCoreBlockPower class

    end

end % end constructors


%% Methods
methods
    function state = iterate(this)
        state = this;
    end
end % end methods

end % end classdef

