classdef BlockPowerAlgoState < handle
%BLOCKPOWERALGOSTATE Current iteration of a block power algorithm
%
%   Contains the current state of a block power iteration algorithm.
%   The state is defined by the vector, and eventually the residual from
%   previous state (equal to NaN in case of first iteration).
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
    % the block vector representing the current solution
    vector;

    % the residual from previous state. NaN for first iteration.
    residual = NaN;

end % end properties


%% Constructor
methods
    function this = BlockPowerAlgoState(init, varargin)
        % Constructor for BlockPowerAlgo class
        %
        %   usage:
        %   STATE = BlockPowerAlgoState(V);
        %
        
        % copy constructor
        if isa(init, 'BlockPowerAlgoState')
            this.vector = init.vector;
            this.residual = init.residual;
            return;
        end
        
        % initialisation constructor
        if ~isa(init, 'AbstractBlockMatrix')
            error('Requires an instance of BlockMatrix as first input');
        end
        this.vector = init;
        
        % eventuallyu copy residual
        if ~isempty(varargin) && isnumeric(varargin{1}) && isscalar(varargin{1})
            this.residual = varargin{1};
        end
    end

end % end constructors

end % end classdef

