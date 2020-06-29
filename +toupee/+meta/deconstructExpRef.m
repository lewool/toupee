function [subject, expDate, expNum, expSpecs] = deconstructExpRef(expRef)
% Deconstructs a Rigbox experiment reference ("exp ref") char array
%
% An expRef is a char array of form '<expDate>_<expNum>_<subject>'.
%
% Inputs:
% -------
% expRef : char array
%  The exp ref, taking the form of '<expDate>_<expNum>_<subject>'.
%
%
% Outputs:
% --------
% subject : char array
%   The experiment's subject.
% expDate : char array
%   The experiment's date.
% expNum : int scalar
%   The experiment's session number.
% expSpecs : cell array
%   Contains `subject`, `expDate`, and `expNum` as cells.
%
% Examples:
% ---------
% 1) Deconstruct exp ref '2020-02-03_1_LEW031'.
%   [subject, expDate, expNum, expSpecs] =...
%       toupee.meta.deconstructExpRef('2020-02-03_1_LEW031');
%

import toupee.misc.extractCellVals

% Ensure `expRef` `isExpRef`
if ~isExpRef(expRef)
    error('toupee:meta:deconstructExpRef:badInput',...
          'The "expRef" input arg is not a valid expRef');
end
expSpecs = strsplit(expRef, '_');
[expDate, expNum, subject] = extractCellVals(expSpecs);
expNum = str2double(expNum);

end