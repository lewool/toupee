function t = cellifyTableVarsNested(t)
% Turns datatypes of a table and its nested tables to cells.
%
% Turns all non-table and non-struct datatypes of a table and its nested
% tables to cells.
%
%
% Inputs:
% -------
% t : table
%   A table, which may contain nested tables.
%
%
% Outputs:
% --------
% t : table
%   The updated table, whose non-table and non-struct datatypes are now all
%   cells.
%
%
% Examples:
% ---------
% 1)
%

% Import self for recursive call.
import toupee.misc.cellifyTableVarsNested
import toupee.misc.iif

% Procedure:
%
%  For each column, check if its a table. 
%   If it is, make a recursive call.
%   Else check if column is a cell or a struct.
%    If it's not, then convert column datatype to a cell.

% Do some checks on input args.
if ~istable(t)
    error('toupee:misc:cellifyTableVarsNested:badInput',...
          'Input arg must be table');
else
    nCols = width(t);  % number of cols
    %nRows = height(t);  % number of rows
    %for iRow = 1:nRows
        for iCol = 1:nCols  % for each col
            cCol = t.(iCol);  % current col
            %cCol = iif(iscell(cCol), @() cCol{1}, cCol);  % get val if cell
            if istable(cCol)  % if col is table, make recursive call.
                t.(iCol) = cellifyTableVarsNested(t.(iCol));
%                 nestedColsMask =...  % check nested cols
%                     table2array(varfun(@(col) isequal(class(col), 'table'),...
%                                        cCol));
%                 if any(nestedColsMask)  % if nested col is table...
%                     t{iRow, iCol}{1} =...% make recursive call
%                         cellifyTableVarsNested(t{iRow, iCol}{1});
%                 end
            else  % if column isn't table...
                % and its not cell and its not struct
                if ~iscell(cCol) && ~isstruct(cCol)
                    % cellify
                    if size(cCol, 1) > 1
                        t.(iCol) = num2cell(cCol);
                    else
                        t.(iCol) = {cCol};
                    end
                end
            end
        end
    %end
end

end

% t.(iCol) = num2cell(t.(iCol))
% cell(blockTable.(2));
% Error using cell
% Size vector should be a row vector with real elements.
% {blockTable.(5)}
% num2cell(e.('newTrialValues'));