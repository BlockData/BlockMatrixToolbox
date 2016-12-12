classdef BlockPowerAlgo < handle
%BLOCKPOWERALGO Power Iteration algorithms for Block-Matrices
%
%   Class BlockPowerAlgo
%
%   Example
%   BlockPowerAlgo
%
%   See also
%     JacobiBlockPower, GaussBlockPower 

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-01-15,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - BIA-BIBS.


%% Properties
properties
end % end properties


%% Constructor
methods
    function this = BlockPowerAlgo(varargin)
    % Constructor for BlockPowerAlgo class

    end

end % end constructors


%% Methods
methods (Abstract)
    iterate(this)
end % end methods

end % end classdef

