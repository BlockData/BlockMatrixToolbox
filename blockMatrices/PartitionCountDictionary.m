classdef PartitionCountDictionary < handle
%PARTITIONCOUNTDICTIONARY Singleton class for counting integer partition
%
%   The PartitionCountDictionary is an utility class used by the
%   IntegerPartition class. Its aim is to keep the singleton instance of
%   the dictionary used to store number of partitions of each integer.
%
%   Usage
%   dict = getDictionary(PartitionCountDictionary.getInstance());
%   Returns the singleton instance of the dictionary, creating a new one if
%   necessary.
%
%   dict is an instance of containers.Map, with keys being char arrays with
%   format '%d.%d'. The first digit corresponds to the integer to be
%   partitioned, the second digit corresponds to the number of partitions.
%   
%
%   Example
%     % Initialize and fill up the dictionary
%     IntegerPartition.countPartitions(6, 3);
%     % get current dictionary instance
%     dict = getDictionary(PartitionCountDictionary.getInstance());
%     % alternative way:
%     % dict = IntegerPartition.getDictionaryInstance();
%     % ask for the number of partitions of integer "4" in two parts
%     get(dict, '4.2')
%     ans =
%          6
%
%   See also
%     IntegerPartition
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2015-07-17,    using Matlab 8.5.0.197613 (R2015a)
% Copyright 2015 INRA - Cepia Software Platform.

%% Properties
properties
    % the dictionary
    dict;
end % end properties

methods (Access = private)
    function this = PartitionCountDictionary(varargin)
        % private constructor to avoid initialisation
        disp('create dictionary');
        this.dict = containers.Map();
    end
end

methods (Static)
    function inst = getInstance()
        % Return the instance of dictionary storing partition counts
        persistent singleton
        if isempty(singleton)
            singleton = PartitionCountDictionary();
        end
        inst = singleton;
    end
end

methods
    function dict = getDictionary(this)
        dict = this.dict;
    end
end
end