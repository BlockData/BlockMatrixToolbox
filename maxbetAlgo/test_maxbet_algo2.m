% Check convergence of maxbet algo using mutli-block data structures
% 
%   Uses data from Hanafi and Ten Berge (2003)
% 
%   USage:
%   test_maxbet_algo2(input)
%
%   Example
%   test_maxbet_algo2
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-06-14,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.


%% setup data

% main problem matrix data, corresponding to X' * X
data = [...
     45  -20    5    6   16    3 ;...
    -20   77  -20  -25   -8  -21 ;...
      5  -20   74   47   18  -32 ;...
      6  -25   47   54    7  -11 ;...
     16   -8   18    7   21   -7 ;...
      3  -21  -32  -11   -7   70 ;...
    ];

% block dimensions of block matrix
mdims = BlockDimensions({[2 2 2], [2 2 2]});

% create the BlockMatrix data structure
A = BlockMatrix(data, mdims);

% block dimensions of block vector
vdims = BlockDimensions({[2 2 2], 1});


%% First initial vector

disp('Initial vector = ones(6,1)');

U0 = BlockMatrix(ones(6, 1), vdims);

U_hat_1 = blockMaxbet_Au(A, U0);

fu1 = U_hat_1' * A * U_hat_1;
disp(sprintf('f(u) = %f \n', fu1.data)); %#ok<DSPS>

disp(U_hat_1');


%% Second initial vector

disp('Initial vector = [ 1 1 1 -1 -1 -1]''');

U0 = BlockMatrix([ones(3, 1); -ones(3,1)], vdims);

U_hat_2 = blockMaxbet_Au(A, U0);

fu2 = U_hat_2' * A * U_hat_2;
disp(sprintf('f(u) = %f ', fu2.data)); %#ok<DSPS>

disp(U_hat_2');

