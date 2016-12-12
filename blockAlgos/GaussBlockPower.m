classdef GaussBlockPower < BlockPowerAlgo
%GAUSSBLOCKPOWER Gauss algorithm for solving block power algorithms
%
%   Class GaussBlockPower
%
%   Example
%   GaussBlockPower
%
%   See also
%     JacobiBlockPower

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-01-15,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - BIA-BIBS.


%% Properties
properties
    % the blockMatrix representing the problem
    data;
    % the block vector representing the initial solution
    init;
    % the current solution
    solution;
    
end % end properties


%% Constructor
methods
    function this = GaussBlockPower(varargin)
    % Constructor for GaussBlockPower class

    end

end % end constructors


%% Methods
methods
end % end methods

end % end classdef

