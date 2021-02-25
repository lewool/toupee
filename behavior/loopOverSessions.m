%% make a list of mice/experiments you want to analyze

%{
mouseList = {...
    {'LEW031'}...
    {'LEW032'}};
%}

mouseName = {{'LEW031'}};
expList = { ...
    {'2020-02-28',2,[2]}};
 %{   
    {'2020-02-03',1,[1]}...    
    {'2020-02-14',1,[1]}...
    {'2020-02-17',1,[1]}...
    {'2020-02-18',1,[1]}...
    {'2020-02-25',1,[1]}...
    {'2020-02-26',1,[1]}...
    {'2020-02-28',2,[2]}...
    };

mouseName = {{'LEW032'}};
expList = { ...
    {'2020-02-03',1,[1]}...
    {'2020-02-13',1,[1]}...
    {'2020-02-14',3,[3]}...
    {'2020-02-17',2,[2]}...
    {'2020-02-18',1,[1]}...
    };
 %}   


%% load all the experiments into expInfo
% process the usual data
% the script knows to loop over all the experiments you listed above
% this will take a while but the command line will print progress
expInfo = initExpInfo(mouseName,expList);
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

%%
%indexes trials by early and late 1st move and makes vectors for all
%Facemap ROIs
%plot graphs of pupil,whisking, paw and neural activity divided by early vs
%late trials
[earlySessionsWhisk,lateSessionsWhisk] = earlyVsLate(expInfo,behavioralData,neuralData,eyeData);

%% Pre-stimulus whisking analysis: compute linear fit and plot whisking +fit
%plot quantification of all early vs late trial pre-stim mean whisk & whisk slope
%needs to first run earlyVsLate function
prestimWhiskAnalysis(eyeData,earlySessionsWhisk,lateSessionsWhisk)
 
%% Rasters for all trials continously sorted (individual cells/ROIs)
%placing all trials in whichTrials  
whichTrials{1}=[];
for t = 1:length(eyeData.eta.alignedFace{1}(:,1,1))
    whichTrials{1}(end+1)=t;
end

%choosing cells correlating to whisking to plot only those
interp_whiskTrace = interp1(eyeData.timeAligned, eyeData.proc.face{1, 2}.motion, neuralData(1).respTimes); 
%remove the NaNs
interp_whiskTrace = interp_whiskTrace(~isnan(interp_whiskTrace));
whiskCells=[];
%start correlating from the time wihtout NaNss in the interpolated
%whisktrace
for c=1:length(neuralData.cellResps(49:end,:))
    [r, p] = corrcoef(interp_whiskTrace, neuralData.cellResps(49:end,c));
    if p(1,2) < 0.05
        whiskCells(end+1)= c;
    end
end

%plot raster of neural activity sorted why pre-stim whisking (200-0 ms)



%% Rasters for grouping trials (individual cells/ROIs)

%grouping trials by high vs low whisking (whisk-sorted) 
[relativeTimes,sortIdxWhisk] = sortTrialByWhisk(whichTrials,eyeData,et,wm);
transpose(sortIdxWhisk);
%high shisking
whiskTrials{1} = (sortIdxWhisk(end-99:end));
%low
whiskTrials{2} = (sortIdxWhisk(1:100,1));

%filter cells that have significant activity during whiskTrials but not during nowhiskTrials 
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);
whiskCells=[];
for icell = 1:length(neuralData.eta.alignedResps{1,3})
    [p,h]=ranksum(baselineResps(whiskTrials,icell),baselineResps(nowhiskTrials,icell));
    if h==1
        whiskCells(end+1)= icell;
    end
end
%choose your Raster plotting parameters:
%whichTrials=
whichCells = whiskCells;
whichROIs = 1:4;
whichSort = 'byWhisk'; %options: 'byWhisk' or 'byEventTime(HAS ERROR THAT NEEDS FIX!)
k=1;
%neural data: all trials in 1 group, 1 psth line
plotWhiskRasters(expInfo, behavioralData, neuralData, whiskCells, whichTrials,whichSort,k)
%Facemap outputs: all trials in 1 group, 1 psth line
plotFacemapRasters(expInfo, behavioralData, eyeData, 1:5, whichTrials,whichSort,k)

%% session raster for population neural activity compared to whisk-trace
%whole session 


%shorter time of sessison 


