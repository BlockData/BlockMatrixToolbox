classdef BlockTable < handle
%BLOCKTABLE Data table of Block-partitioned data
%
%   (Uses block partition only for columns)
%
%   Example
%   BlockTable
%
%   See also
%

% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2016-01-18,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2016 INRA - BIA-BIBS.


%% Properties
properties
    % inner data of the table, stored in a Nr-by-Nc BlockMatrix.
    % The block dimension of the block-matrix will be use to initialize the
    % partition of row names and column names.
    data;
    
    % name of columns, stored in a 1-by-Nc cell array of strings
    colNames;
    
    % name of rows, stored in a Nr-by-1 cell array of strings
    rowNames;
    
    % the names of block-columns, stored on a 1-by-NBC cell array of strings
    % where NBC is the number of blocks in horizontal direction
    blockColNames;
    
    % the names of block-rows, stored on a NBR-by-1 cell array of strings
    % where NBR is the number of blocks in vertical direction
    blockRowNames;
    
end % end properties


%% Constructor
methods
    function this = BlockTable(varargin)
        % Constructor for BlockTable class
        
        if isa(varargin{1}, 'BlockTable')
            % copy constructor
            bt = varargin{1};
            this.data = bt.data;
            this.rowNames   = bt.rowNames;
            this.colNames   = bt.colNames;
            
            varargin(1) = [];

        elseif isa(varargin{1}, 'BlockMatrix')
            % Default constructor, using BlockMatrix as first input,
            % possibly followed by other arguments
            this.data = varargin{1};
            varargin(1) = [];
           
        end
        
        % check if column names were specified
        if ~isempty(varargin)
            if iscell(varargin{1})
                if length(varargin{1}) ~= size(this.data,2)
                    error('Column names have %d elements, whereas table has %d columns', ...
                        length(varargin{1}), size(this.data,2));
                end
                this.colNames = varargin{1};
                varargin(1) = [];
            end
        end
        
        % check if row names were specified
        if ~isempty(varargin)
            if iscell(varargin{1})
                if length(varargin{1}) ~= size(this.data,1)
                    error('Number of row names does not match row number');
                end
                this.rowNames = varargin{1};
                varargin(1) = [];
            end
        end

        % ---------
        % parse additional arguments set using parameter name-value pairs
        
        while length(varargin) > 1
            % get parameter name and value
            param = lower(varargin{1});
            value = varargin{2};
            
            % switch
            if strcmp(param, 'rownames')
                if length(value) ~= size(this.data,1)
                     error('Number of row names does not match row number');
                end
                this.rowNames = value;
                
            elseif strcmp(param, 'colnames')
                if length(value) ~= size(this.data,2)
                     error('Number of column names does not match column number');
                end
                this.colNames = value;

            elseif strcmp(param, 'blocknames')
                if length(value) ~= blockSize(this.data, 2)
                     error('Number of block names does not match block-column number');
                end
                this.blockColNames = value;

            
            else
                error('BlockTable:BlockTable', ...
                    ['Unknown parameter name: ' varargin{1}]);
            end
            
            varargin(1:2) = [];
        end
        
        % ---------
        % create default values for other fields if they're not initialised
        
        % size if the data table
        nr = size(this.data, 1);
        nc = size(this.data, 2);
        
        if isempty(this.rowNames) && nr > 0
            this.rowNames = strtrim(cellstr(num2str((1:nr)')));
        end
        if isempty(this.colNames) && nc > 0
            this.colNames = strtrim(cellstr(num2str((1:nc)')))';
        end
        
        % block-size of data block-matrix
        nbr = blockSize(this.data, 1);
        nbc = blockSize(this.data, 2);
        
        if isempty(this.blockRowNames) && nbr > 0
            this.blockRowNames = strtrim(cellstr(num2str((1:nbr)')));
        end
        if isempty(this.blockColNames) && nbc > 0
            this.blockColNames = strtrim(cellstr(num2str((1:nbc)')))';
        end
        
    end

end % end constructors


%% Display methods
methods
    function disp(this, varargin)
        
        %% Initialisations
        
        % loose format: display more empty lines
        isLoose = strcmp(get(0, 'FormatSpacing'), 'loose');
        
        % isLong = ~isempty(strfind(get(0,'Format'),'long'));
        % dblDigits = 5 + 10*isLong; % 5 or 15
        % snglDigits = 5 + 2*isLong; % 5 or 7
%         maxWidth = get(0, 'CommandWindowSize');
%         maxWidth = maxWidth(1);
        
        % get BlockMatrix total size
        dim = size(this.data);
        nRows = dim(1);
        nCols = dim(2);
        
        % process case of empty table
        if nRows == 0 || nCols == 0
            disp('empty BlockTable');
            return;
        end
        
        % Display information on block sizes
        className = class(this);
        disp(sprintf('%s object with %d rows and %d columns', className, nRows, nCols)); %#ok<DSPS>

        % get BlockMatrix total size
        dim = blockSize(this.data);
        nBlockCols = dim(2);
        
        if isLoose
            fprintf('\n');
        end
        
        % extract the block partition along columns
        colParts = blockDimensions(this.data, 2);
        
        % display constants
        nSpacesBetweenCols = 2;
        nSpacesBetweenBlocks = 4;
        
        % data access shortcuts
        data0 = this.data.data;
        
        % list of formats for the content of each column
        colFormats = cellstr(repmat('%2d', nCols, 1))';
        rowFormat = '%s';
        for iCol = 1:nCols
            rowFormat = [rowFormat ' ' colFormats{iCol}]; %#ok<AGROW>
        end

        % determine width of each data column
        colWidths = zeros(1, nCols);
        dataColWidths = zeros(1, nCols);
        for iCol = 1:nCols
            % determine data column width by trying to format data and
            % measuring resulting length
            w = 0; 
            fmt = colFormats{iCol};
            for i = 1:size(data0, 1)
                w = max(w, length(sprintf(fmt, data0(i, iCol)))); 
            end
            dataColWidths(iCol) = w;
            
            % determine maximum width of data and of column names
            colWidths(iCol) = max(w, length(this.colNames{iCol}));
        end

        rowNamesWidth = 8;
        rowNamesBlanks = repmat(' ', 1, rowNamesWidth);
        
        % display names of each block column
        iCol = 0; % index of column
        fprintf(rowNamesBlanks);
        for iColBlock = 1:nBlockCols
            % display each column of the current block
            nCols2 = colParts(iColBlock);
            
            % first column of the block
            iCol = iCol + 1;
            nChars = colWidths(iCol);
            
            % remaining columns of the block
            for iBlockCol = 2:nCols2
                iCol = iCol + 1;
                nChars = nChars + nSpacesBetweenCols + colWidths(iCol);
            end
            
            fprintf(repmat(' ', 1, nSpacesBetweenBlocks));
            fmt = ['%-' num2str(nChars) 's'];
            fprintf(fmt, this.blockColNames{iColBlock});
        end
        fprintf('\n');
        
        
        % display small lines to mark block columns
        iCol = 0; % index of column
        fprintf(rowNamesBlanks);
        for iColBlock = 1:nBlockCols
            % display each column of the current block
            nCols2 = colParts(iColBlock);
            % first column of the block
            iCol = iCol + 1;
            fprintf(repmat(' ', 1, nSpacesBetweenBlocks));
            fprintf(repmat('_', 1, colWidths(iCol)));
            % remaining columns of the block
            for iBlockCol = 2:nCols2
                iCol = iCol + 1;
                fprintf(repmat('_', 1, nSpacesBetweenCols + colWidths(iCol)));
            end
        end
        fprintf('\n');
        
        
        % display column names
        iCol = 0; % index of column
        fprintf(rowNamesBlanks);
        for iColBlock = 1:nBlockCols
            % display each column of the current block

            % spacing before first column
            nCols2 = colParts(iColBlock);
            fprintf(repmat(' ', 1, nSpacesBetweenBlocks));
            
            % first column of the block
            iCol = iCol + 1;
            fmt = ['%' num2str(colWidths(iCol)) 's'];
            fprintf(fmt, this.colNames{iCol});
            
            % remaining columns of the block
            for iBlockCol = 2:nCols2
                iCol = iCol + 1;
                fprintf(repmat(' ', 1, nSpacesBetweenCols));
                fmt = ['%' num2str(colWidths(iCol)) 's'];
                fprintf(fmt, this.colNames{iCol});
            end
        end
        fprintf('\n');

        % display column names
        iCol = 0; % index of column
        fprintf(rowNamesBlanks);
        for iColBlock = 1:nBlockCols
            % display each column of the current block

            % spacing before first column
            nCols2 = colParts(iColBlock);
            fprintf(repmat(' ', 1, nSpacesBetweenBlocks));
            
            % first column of the block
            iCol = iCol + 1;
            fprintf(repmat('_', 1, colWidths(iCol)));
            
            % remaining columns of the block
            for iBlockCol = 2:nCols2
                iCol = iCol + 1;
                fprintf(repmat(' ', 1, nSpacesBetweenCols));
                fprintf(repmat('_', 1, colWidths(iCol)));
            end
        end
        fprintf('\n');

        
        % iterate over data rows
        for iRow = 1:nRows
            % display row name
            fmt = ['%' num2str(rowNamesWidth) 's'];
            fprintf(fmt, this.rowNames{iRow});
            
            % iterate over block-columns
            iCol = 0; % index of column
            for iColBlock = 1:nBlockCols
                % display each column of the current block
                
                % first column of the block
                iCol = iCol + 1;

                % spacing before first column
                nCols2 = colParts(iColBlock);
                fprintf(repmat(' ', 1, nSpacesBetweenBlocks));
                
                % display value of first column
                value = data0(iRow, iCol);
                str = sprintf(colFormats{iCol}, value);
                fmt = ['%-' num2str(colWidths(iCol)) 's'];
                fprintf(fmt, str);
                
                % remaining columns of the block
                for iBlockCol = 2:nCols2
                    iCol = iCol + 1;
                    
                    % spacing with previous column
                    fprintf(repmat(' ', 1, nSpacesBetweenCols));
                    
                    % display value of first column
                    value = data0(iRow, iCol);
                    str = sprintf(colFormats{iCol}, value);
                    fmt = ['%-' num2str(colWidths(iCol)) 's'];
                    fprintf(fmt, str);
                end
            end
            fprintf('\n');
            
        end
        
        % final formatting
        if isLoose
            fprintf('\n');
        end

    end
end % end display methods
    
end % end classdef

