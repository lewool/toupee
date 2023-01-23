function [kernelAnalysis] = loadKernelAnalysis(mouseName, expDate, expNum)

sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));

m = load(fullfile(sessDir,'kernelAnalysis.mat'));

behavioralData = m.behavioralData;
expInfo = m.expInfo;
neuralData = m.neuralData;
