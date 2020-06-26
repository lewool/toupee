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
    for fnameIdx = 1:numel(fnames)
        curFld = s.(fnames{fnameIdx});
        if isequal(class(curFld), 'struct')
            curFldFnames = fieldnames(curFld);
            curFldFnamesStructMask =...
                cellfun(@(fname) isequal(class(curFld.(fname)), 'struct'),...
                        curFldFnames);
            if any(curFldFnamesStructMask)
                s.(fnames{fnameIdx}) = struct2tableNested(curFld);
            else
                s.(fnames{fnameIdx}) = struct2table(curFld);
            end
        end
    end
    t = struct2table(s);
end

end