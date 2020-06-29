function t = struct2tableNested(s)
% Turns a struct and all its nested structs into tables
% 
%
% Inputs:
% -------
% s : struct
%   A struct, which may contain nested structs.
%
%
% Outputs:
% --------
% t : table
%   A table, which may contain nested tables.
%
%
% Examples:
% ---------
% 1) Convert a triple-nested struct into a triple-nested table:
%   s.a = 1;
%   s.b = struct();
%   s.c = struct();
%   s.b.e = 1;
%   s.b.f = struct();
%   s.b.g = struct();
%   s.b.f.h = 1;
%   s.b.g.i = 1;
%   s.c.j = 1;
%   s.c.k = struct();
%   s.c.l = struct();
%   s.c.k.m = 1;
%   s.c.l.n = 1;
%   t = toupee.misc.struct2tableNested(s);
%


% Import self for recursive call.
import toupee.misc.struct2tableNested

% Procedure:
% For each fieldname, check if its a struct. 
%   If it is, check its fieldnames. 
%       If any of those nested fieldnames is a struct, make recursive call.
%       Else convert fieldname to table and add it to struct.
%
% Convert struct to table.

% Do some checks on input args.
if ~isstruct(s)
    error('toupee:misc:struct2tableNested:badInput',...
          'Input arg must be struct');
else
    fnames = fieldnames(s);
    for iFname = 1 : numel(fnames)  % for each fieldname
        curFld = s.(fnames{iFname});  % get the field
        if isequal(class(curFld), 'struct')  % check if its a struct
            % See if any of the field's fields are structs
            curFldFnames = fieldnames(curFld);
            curFlds = cellfun(@(fname) curFld.(fname), curFldFnames,...
                              'uni', 0);
            curFldFnamesStructMask =...
                cellfun(@(fld) isequal(class(fld), 'struct'), curFlds);
            % if field's fields are structs, make recursive call
            if any(curFldFnamesStructMask)
                s.(fnames{iFname}) = struct2tableNested(curFld);
            else  % else convert the field to a table
                % for each struct in the field's struct array
                for iS = 1:numel(s)  
                    % If struct fields have different lengths, create
                    % single-row tables
                    try
                        s(iS).(fnames{iFname}) = struct2table(curFld);
                    catch ex
                        if strcmp(ex.identifier,...
                                  ['MATLAB:struct2table:UnequalField'...
                                   'Lengths'])
                            s(iS).(fnames{iFname}) =...
                                struct2table(curFld, 'AsArray', 1);
                        else
                            rethrow(ex);
                        end
                    end
                end
            end
        end
    end
    % If struct fields have different lengths, create single-row tables
    try
        t = struct2table(s);
    catch ex
        if strcmp(ex.identifier, 'MATLAB:struct2table:UnequalFieldLengths')
            t = struct2table(s, 'AsArray', 1);
        else
            rethrow(ex);
        end
    end
end

end

