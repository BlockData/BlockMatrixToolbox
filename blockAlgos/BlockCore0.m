classdef BlockCore0 < BlockCore
%BLOCKCORE0 first demo of Blockcore implementation
%
%   Class BlockCore0
%
%   Example
%   BlockCore0
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
    data;
end % end properties


%% Constructor
methods
    function this = BlockCore0(matrix)
    % Constructor for BlockCore0 class
        this.data = matrix;
    end

end % end constructors


%% Methods
methods
    function res = getMatrix(this, vector)
        res = this.data;
    end
end % end methods

end % end classdef

