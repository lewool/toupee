function [expRef, expLog] = constructExpRef(subject, expDate, expNum)
% Constructs a Rigbox experiment reference ("exp ref") char array
%
% An expRef is a char array of form '<expDate>_<expNum>_<subject>'.
%
%
% Inputs:
% -------
% subject : char array
%   The experiment's subject.
% expDate : char array
%   The experiment's date. Must be in form 'yyyy-mm-dd'.
% expNum : int scalar
%   The experiment's session number.
%
%
% Outputs:
% --------
% expRef : char array
%  The exp ref, taking the form of '<date>_<session>_<subject>'.
% expLog : cell array
%  The subject, expDate, and expNum as separate entries in a cell.
%
%
% Examples:
% ---------
% 1) Create exp ref for 'LEW031', '2020-02-03', 1.
%   [expRef, expLog] =...
%       toupee.meta.constructExpRef('LEW031', '2020-02-03', 1);
%

% Import all other functions in this subpackage.
import toupee.meta.*
% Make expRef.
expRef = sprintf('%s_%i_%s', expDate, expNum, subject);
% Check to make sure `expRef` is a valid expRef.
if ~isExpRef(expRef)
    warning('toupee:meta:constructExpRef:badExpRef',...
            ['The constructed expRef, "%s", is not a valid expRef. The ',...
             'input args should be a string, a string in a "yyyy-mm-dd" ',...
             'format, and an int scalar. The input args were: "%s" "%s" ',...
             '"%i'], expRef, subject, expDate, expNum); 
end
expLog = {subject, expDate, num2str(expNum)};

end