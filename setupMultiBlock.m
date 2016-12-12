function setupMultiBlock(varargin)
%SETUPBLOCKBERRY Setup paths required to run the BlockBerry library
%
%   usage:
%   setupBlockBerry
%   
%   Get start by typing 'help multiBlock'
%
 
% ------
% Author: David Legland
% e-mail: david.legland@nantes.inra.fr
% Created: 2015-04-30,    using Matlab 8.3.0.532 (R2014a)
% Copyright 2015 INRA - Cepia Software Platform.

% extract parent path of the program
fileName = mfilename('fullpath');
mainDir = fileparts(fileName);

disp('Installing BlockMatrix Toolbox...');

% add libraries
addpath(fullfile(mainDir, 'blockTensors'));
addpath(fullfile(mainDir, 'blockMatrices'));
addpath(fullfile(mainDir, 'blockAlgos'));
addpath(fullfile(mainDir, 'blockTables'));

disp('BlockMatrix Toolbox installed!');

