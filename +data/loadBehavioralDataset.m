function [behavioralData, expInfo] = loadBehavioralDataset(mouseName, expDate, expNum)

sessDir = fullfile('G:\Analysis\Behavior',mouseName,expDate,num2str(expNum));

m = load(fullfile(sessDir,'dataset.mat'));

behavioralData = m.behavioralData;
expInfo = m.expInfo;