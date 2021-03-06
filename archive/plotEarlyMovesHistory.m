% mouseList = ...
%     {{'LEW005'},...
%     {'LEW006'},...
%     {'LEW006'},...
%     {'LEW006'},...
%     {'LEW013'},...
%     {'LEW013'},...
%     {'LEW015'},...
%     {'LEW015'},...
%     {'LEW015'}};
% 
% expList = ...
%     {{'2018-06-10',2,[2 3]},...
%     {'2018-06-13',2,[2 3]},...
%     {'2018-06-14',1,[1 2]},...
%     {'2018-06-15',1,[1 2]},...
%     {'2019-03-26',1,1},...
%     {'2019-03-27',1,1},...
%     {'2019-03-19',1,1},...
%     {'2019-03-21',1,1},...
%     {'2019-04-12',1,1}};

mouseList = ...
    {{'LEW015'}};

expList = ...
    {{'2019-03-18',1,1},...
    {'2019-03-19',1,1},...
    {'2019-03-21',1,[2 3]},...
    {'2019-03-22',1,[1 2]},...
    {'2019-03-25',1,[1 2]},...
    {'2019-03-26',2,1},...
    {'2019-03-27',1,1},...
    {'2019-03-29',1,1},...
    {'2019-04-01',1,1},...
    {'2019-04-02',1,1},...
    {'2019-04-03',1,1},...
    {'2019-04-04',2,1},...
    {'2019-04-09',1,1},...
    {'2019-04-10',1,1},...    
    {'2019-04-12',1,1}};

%%    

expInfo = initExpInfo(mouseList,expList);

%% load data
for i = 1:length(expList)
    expInfo(i) = data.loadExpData(expInfo(i));
    [eventTimes{i}, wheelTrajectories{i}] = getEventTimes(expInfo(i), {'stimulusOnTimes' 'interactiveOnTimes' 'stimulusOffTimes'});
end

%% 
firstMoveTimes = [];
interactiveOnTimes = [];
earlyMoves = [];

for i = 1:length(expList)
    firstMoveTimes = cat(2,firstMoveTimes, eventTimes{i}(7).daqTime);
    interactiveOnTimes = cat(2,interactiveOnTimes, eventTimes{i}(2).daqTime);
    sessionEnds(i) = length(firstMoveTimes);
end
earlyMoves = firstMoveTimes - interactiveOnTimes <= 0;

%%

plot(movmean(earlyMoves,100));
hold on;
for e = 1:length(sessionEnds)-1
line([sessionEnds(e) sessionEnds(e)],[0 1]);
end
