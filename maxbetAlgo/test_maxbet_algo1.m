% Check convergence of maxbet algo using specific function maxbet_Au
%
%   Uses data from Hanafi and Ten Berge (2003)
%
% Usage:
%   test_maxbet_algo1
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-06-13,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - Cepia Software Platform.


%% setup data

% main problem matrix, corresponding to X' * X
A = [...
     45  -20    5    6   16    3 ;...
    -20   77  -20  -25   -8  -21 ;...
      5  -20   74   47   18  -32 ;...
      6  -25   47   54    7  -11 ;...
     16   -8   18    7   21   -7 ;...
      3  -21  -32  -11   -7   70 ;...
    ];

% block size along each dimension
vdims = [2 2 2];


%% First initial vector

disp('initial vector = ones(6,1)');

u0 = ones(6, 1);

u_hat_1 = maxbet_Au(A, u0, vdims);

disp(u_hat_1');

fu1 = u_hat_1' * A * u_hat_1;
disp(sprintf('f(u) = %f \n', fu1)); %#ok<DSPS>


%% Second initial vector

disp('initial vector = [ 1 1 1 -1 -1 -1]''');

u0 = [ones(3, 1); -ones(3,1)];

u_hat_2 = maxbet_Au(A, u0, vdims);

disp(u_hat_2');

fu2 = u_hat_2' * A * u_hat_2;
disp(sprintf('f(u) = %f ', fu2)); %#ok<DSPS>
