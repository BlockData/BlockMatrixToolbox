function coreModel = core_methods(data, method, maxdim, init, CORE, tol)
%CORE_METHODS implement "core methods" for block matrices
%
%   COREMODEL = core_methods(DATA, METHOD, MAXDIM, INIT, CORE, TOL)
%
%   Input arguments:
%   DATA:       the input data table, as a multiblock matrix
%   METHOD:     the name of the method to use. Typically: 'maxbet'
%   MAXDIM:     the maximal number of dimension to use for factorisation (?)
%   INIT:       initial value of the solution (?)
%   CORE:       the type of deflation used (can be on global scores, on
%               block loadings, on blck scores...). One of 'CORE1',
%               'CORE2', 'CORE3'.
%   TOL:        the tolerance used for computing convergence
%
%   Output argument:
%   CORE:       an instance of CoreModel, representing a decomposition of
%               the input block matrix
%
%   Example
%   core_methods
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-04-14,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - Cepia Software Platform.

%% Initialisation

% add test in the data structure
data0 = data;

% compute the largest possible value of maxdim
% TODO: implement function
m = optimal_maxdim(data0, CORE);

% ensure maxdim is lower or equal to the largest possible maxdim value
maxdim = min(m, maxdim);

% number of blocks in the BlockMatrix
maxblo = size(data, 2);

% initial value 
tt = init;


%% Main part

switch lower(method)
    case 'maxbet'
        % iteration sur les dimensions de reduction
        for dim = 1:maxdim
            % concatenation des u_i
           q = maxbet_procedure1(data0, tt, tol);
            % si on veut choisir le type de maxbet, utiliser la ligne
            % suivante:
            % q = maxbet_procedures(data0, tt, tol, 1);
            
            % calcule les parametres pour la reduction de dimension
            % TODO: implement function
            [a, b] = deflation_parameters(data0, q, CORE);
            
            % TODO: allocate memory for A and B
            
            switch CORE
                case 'CORE1'
                    % deflate on global scores
                    A{1}(:, dim) = a{1};
                    B{1}(:, dim) = b{1};
                    
                case 'CORE2'
                    % deflate on block loadings
                    for blo = 1:maxblo
                        A{blo}(:, dim) = a{blo};
                        B{1}(:, dim) = b{1};
                    end
                    
                case 'CORE3'
                    % deflate on block scores
                    for blo = 1:maxblo
                        A{blo}(:, dim) = a{blo};
                        B{blo}(:, dim) = b{blo};
                    end

                otherwise
                    error(['Unrecognized method name: ' method]);
            end
            
            % Applique l'algorithme de reduction de dimension
            % TODO: implement function
            % TODO: maybe return result in varaiable 'data' ?
            data0 = WedderburnReduction(data0, a, b, CORE);
        end
        
        % calcule la factorisation
        % TODO: implement function
        coreModel = procedure_core_by_block(data, A, B, CORE);
        
    case lower('ACOM')
        error('method ACOM is not yet implemented');
        
    otherwise
        error(['Unrecognized method name: ' method]);
end

