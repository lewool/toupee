function [behavioralData, expInfo, neuralData] = loadDataset(mouseName, expDate, expNum)

sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));

m = load(fullfile(sessDir,'dataset.mat'));

behavioralData = m.behavioralData;
expInfo = m.expInfo;
neuralData = m.neuralData;
