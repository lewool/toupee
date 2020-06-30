% Advantages of using tables over struct arrays: Can more easily access and
% manipulate data -> more convenient indexing. E.g. can more easily combine
% data across sessions and/or across data vars.

% sessions = rows; data vars = columns;

% can index with names (char arrays) or numbers (int arrays_
% indexing in with `()`, `{}`, or `.`, -> get subtables, cell arrays, or
% underlying values within data vars: 
% expInfo(rows, cols) -> subtable
% expInfo.(col) -> cell array of full col data
% expInfo{rows, cols} -> concatenated cell array of rows, cols data
% expInfo.(col){row} -> underlying data of row, col
% expInfo{row, col}{1} -> underlying data of row, col

% get single data var from single session

% get multiple data vars from single session

% get single data var from multiple sessions

% get multiple data vars from multiple sessions

% add and delete sessions and data vars

% assign into an existing or new data var

% sort/reorder sessions and data vars (sessions: by date or by subject. 
% data vars: alphabetically)

% convert text in data var -> categorical data var (example add expDef as
% data var to `expInfo`)

% combine and split tables

% create custom table properties

% computations & plots with tables