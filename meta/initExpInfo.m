function expInfo = initExpInfo(varargin)

% SINGLE MOUSE (MATCHED EXPS)
% varargin{1} = 'LEW007';

% SINGLE EXP
% varargin{1} = {{'LEW007'}};
% varargin{2} = {'2018-10-09',1,[1]}; 

% SINGLE MOUSE, SINGLE EXP
% varargin{1} = {{'LEW008'}};
% varargin{2} = ...
%     {{'2019-01-29',1,[1]},...
%     {'2019-02-07',1,[1]},...
%     {'2019-02-12',1,[1]}}; 

% MULTI MOUSE, MULTI EXP
% varargin{1} = ...
%     {{'LEW005'},...
%     {'LEW015'},...
%     {'LEW015'}};
% varargin{2} = ...
%     {{'2018-06-10',2,[2 3]},...
%     {'2019-03-21',1,1},...
%     {'2019-04-12',1,1}};

%% check arguments

if nargin == 1 %matched-cell experiments
    
    mouseName = varargin{1};
    paths = data.dataPaths();
    procDir = paths.local{1}{1};
    if ~exist(char(fullfile(procDir,mouseName,'matchedCells')),'dir')
        error('No match folder found for this mouse')
    end

    load(fullfile(procDir,mouseName,'matchedCells','plane1','matchFile.mat'));
    numSessions = length(roiMatchData.allRois);
    
    %get expList straight from matchCell data
    for m = 1:numSessions
        fileParts = regexp(roiMatchData.allRois{m},'\','split');
        expList{m}{1} = fileParts{5};
        expList{m}{2} = str2num(fileParts{6});
        expList{m}{3} = str2num(fileParts{6});
    end

    mouseList = {{mouseName}};
    matchTag = 1;
elseif nargin == 2
    
    mouseList = varargin{1};
    expList = varargin{2};
    matchTag = 0;

elseif nargin > 2
    error('Inputs should be a mouse name or a mouseList/expList')
end
%% assemble expInfo struct

if length(mouseList) == 1
    mouseList = repelem(mouseList,length(expList));
end

    for s = 1:length(expList)
        expInfo(s) = struct(...
            'mouseName',mouseList{s}{:},...
            'expDate',expList{s}{1},...
            'expNum',expList{s}{2},...
            'expSeries',expList{s}{3},...
            'cellMatched',matchTag,...
            'block',[],...
            'Timeline',[],...
            'numPlanes',[],...
            'numChannels',[]);
    end
end