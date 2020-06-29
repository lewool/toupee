function tf = isExpRef(a)
% Checks if input arg is a Rigbox expRef.
%
% An expRef is a char array of form '<expDate>_<expNum>_<subject>'.
%
%
% Inputs:
% -------
% a : any datatype
%
%
% Outputs:
% --------
% tf : logical
%   True/False. `true` if `a` is an expRef, `false` otherwise.
%
%
% Examples:
% ---------
% 1) Check if the following char array is an actual expRef:
%   toupee.meta.isExpRef('2020-02-03_1_LEW031')
%
% See Also:
% ---------
% toupee.meta.constructExpRef
% toupee.meta.deconstructExpRef
%

% Must be char array, must have two underscores, must start with 10 element
% date char array, and must have an integer value for expNum.

% Initialize `tf` as false; set to true if meets all conditions.
tf = false;
% If there is an error trying to meet any condition, ignore error.
try
    if...
    ischar(a) ...
    && numel(find(a == '_')) == 2 ...
    && find(a == '_', 1) == 11 ...
    && numel(num2str(str2num(a(1:10)))) == 4 ...
    && mod(str2num(a(12:find(a == '_', 1, 'last') - 1)), 1) == 0
        tf = true;
    end  %#ok<*ST2NM> 
catch
end

end