classdef BlockCore < handle
%BLOCKCORE Interface for kernel of BlockCore Power algorithm
%
%   Class BlockCore
%
%   Example
%   BlockCore
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
end % end properties


%% Constructor
methods
    function this = BlockCore(varargin)
    % Constructor for BlockCore class

    end

end % end constructors


%% Methods
methods (Abstract)
    getMatrix(this, vector)
        
end % end methods

end % end classdef

