%TEST_MAXBET1 Programme de test pour la fonction 'maxbet_procedure1'
%
%   output = test_maxbet1(input)
%
%   Example
%   test_maxbet1
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-04-14,    using Matlab 8.4.0.150421 (R2014b)
% Copyright 2015 INRA - Cepia Software Platform.

%% Cree une bloc-matrice pour les donnees
% On part sur une bloc-matrice de 3 blocs concatenes en ligne, 
% avec des tailles de bloc de 2, 3 et 2 colonnes, 
% et 4 lignes pour chaque bloc.
% La matrice equivalente a une taille de 4-par-7.

% create block dimensions
mdims = BlockDimensions({4, [2 3 2]});

% create block matrix instance
data = BlockMatrix(rand(4, 7), mdims);

% display the BlockMatrix
disp('Block matrix:');
disp(data);
% displayData(data);


%% Cree une bloc-matrice pour les vecteurs
% On utilise une bloc-matrice de 3 blocs, chaque bloc contenant un vecteur.
% La longueur des vecteurs est de (2,3,2).

% create block dimensions
vdims = BlockDimensions({[2 3 2], 1});

% create block matrix instance
tt = BlockMatrix(rand(7, 1), vdims);

% display the block-vector (transposed)
disp('Transpose of input vector t:');
disp(tt');
% displayData(tt');


%% Run the maxbet procedure

% sample line to run the procedure
q = maxbet_procedure1(data, tt, .01);

% display the result (transposed)
disp('Transpose of result q:');
disp(q');
% displayData(q')
