function m = mergeObj(objs)
% Merges multiple objects into one
%
% Keeps the values for the fields of the objects specified by their reverse
% order into the function. E.g. if we have two objects `a` and `b` that
% each have a field `c`, then for `m = mergeObj(a, b)`, `m.c` will take the
% value of `b.c` instead of `a.c`.
%
%
% Inputs:
% -------
% objs : cell array of structs OR tables
%   Objects to merge.
%
%
% Outputs:
% --------
% m : struct OR table
%   Single, merged object.
%
%
% Examples:
% ---------
% 1) Merge three nested structs:
%   s1.a = 1;
%   s1.b = struct();
%   s1.c = struct();
%   s1.b.d = 1;
%   s1.c.e = struct();
%   s1.c.e.f = 1;
%   s2.a = 2;
%   s2.b = struct();
%   s2.b.e = 2;
%   s3.c = struct();
%   s3.c.e = struct();
%   s3.c.e.f = 3;
%   s3.c.e.g = 3;
%   m = toupee.misc.mergeObj({s1, s2, s3});
%

% Import self for recursive call.
import toupee.misc.mergeObj

% Do some checks on input args.
if numel(objs) < 2
    error('toupee:misc:mergeObj:badInput',...
          'Requires at least 2 objects to merge')
elseif numel(objs) > 2
    leftovers = objs(3:end);  % save these for recursive call
else
    leftovers = {};
end
% Merge 2 objects at a time.
a = objs{1};
b = objs{2};
if ~isequal(class(a), class(b))
    error('toupee:misc:mergeObj:badInput',...
          'All input objects must have the same class')
end

switch class(a)
    case 'struct'  % get fieldnames
        m = b;
        fldnmsA = fieldnames(a);
    case 'table'  % get column names
        m = b;
        fldnmsA = a.Properties.VariableNames;
end

% Add each fieldname in `a` to `m`.
for i = 1:numel(fldnmsA)
    m.(fldnmsA{i}) = a.(fldnmsA{i});
end

if ~isempty(leftovers)  % make recursive call
    leftovers{end + 1} = m;
    m = mergeObj(leftovers);
end

end