function expInfo = initExpInfo(details)
% Gets experiment information for given session(s)
% 
% Inputs:
% -------
% details : cell
%   The session(s) details. Use a nested cell for each session, with four
%   elements in each nested cell: subject name, datestr, session, and
%   series.
% 
% Outputs:
% --------
% expInfo : struct
%   A struct array with each element containing fields with information for
%   a particular session. The fields for each struct element are:
%   'mouseName', 'expDate', 'expNum', 'expSeries', 'block', 'Timeline',
%   'numPlanes', 'numChannels'.
%
% Examples:
% ---------
% 1) Return experiment info for a single session from a single subject.
%   details = {{'LEW031', '2020-02-03', 1, 1}};
%   expInfo = toupee.meta.initExpInfo(details);
% 2) Return experiment info for multiple sessions from multiple subjects.
%   details = {{'LEW005', '2018-06-10', 2, [2 3]},... 
%              {'LEW015', '2019-03-21', 1, 1},...
%              {'LEW015', '2019-04-12', 1, 1}};
%   expInfo = toupee.meta.initExpInfo(details);

% Assemble expInfo struct
for s = 1:length(details)
    expInfo(s) = struct(...
        'mouseName', details{s}{1},...
        'expDate', details{s}{2},...
        'expNum', details{s}{3},...
        'expSeries', details{s}{4},...
        'block', [],...
        'Timeline', [],...
        'numPlanes', [],...
        'numChannels', []);  %#ok<AGROW>
end

end