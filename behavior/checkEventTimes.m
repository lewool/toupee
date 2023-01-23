for m = 1:length(mouseList)
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)
    [behavioralData, expInfo] = data.loadBehavioralDataset(mouseName, expDate, expNum);
    [n(m,:), ~] = histcounts(behavioralData.eventTimes(2).daqTime - behavioralData.eventTimes(1).daqTime,linspace(0,2,41));
    mins(m) = min(behavioralData.eventTimes(2).daqTime - behavioralData.eventTimes(1).daqTime);
end