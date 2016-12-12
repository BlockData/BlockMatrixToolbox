classdef CoreModel < handle
%COREMODEL  Represent a factorisation of a block-matrix
%
%   Class CoreModel
%
%   Example
%   CoreModel
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-04-08,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - BIA-BIBS.


%% Properties
properties
    % a m-by-m BlockMatrix 
    U;
    
    % a m-by-n Block Diagonal Matrix
    S;
    
    % a n-by-n BlockMatrix 
    V;
end % end properties


%% Constructor
methods
    function this = CoreModel(varargin)
    % Constructor for CoreModel class

        if nargin == 3
            this.U = varargin{1};
            this.sigma = varargin{2};
            this.Vstar = varargin{3};

        else
            error('Should specify U, S and V arguments');
        end
    end

end % end constructors


%% Methods
methods
end % end methods

end % end classdef

