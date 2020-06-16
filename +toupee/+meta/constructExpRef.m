function [expRef, expLog] = constructExpRef(subject, expDate, expNum)
% Constructs a Rigbox experiment reference ("exp ref") char array
% This exp ref takes the form of '<date>_<session>_<subject>'.
%
% Inputs:
% -------
% subject : char array
%   The experiment's subject.
% expData : char array
%   The experiment's date.
% expNum : int scalar
%   The experiment's session number.
%
% Outputs:
% --------
% expRef : char array
%  The exp ref, taking the form of '<date>_<session>_<subject>'.
% expLog : cell array
%  The subject, expDate, and expNum as separate entries in a cell.
%
% Examples:
% ---------
% 1) Create exp ref for 'LEW031', '2020-02-03', 1.
%   [expRef, expLog] =...
%       toupee.meta.constructExpRef('LEW031', '2020-02-03', 1);

expRef = sprintf('%s_%i_%s', expDate, expNum, subject);
expLog = {subject, expDate, num2str(expNum)};

end