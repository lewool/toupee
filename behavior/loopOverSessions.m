%% make a list of mice/experiments you want to analyze

%{
mouseList = {...
    {'LEW031'}...
    {'LEW032'}};
%}

mouseName = {{'LEW031'}};
expList = {{'2020-02-28',2,[2]}};   

%{ 
    {'2020-02-03',1,[1]}...    
    {'2020-02-14',1,[1]}...
    {'2020-02-17',1,[1]}...
    {'2020-02-18',1,[1]}...
    {'2020-02-25',1,[1]}...
    {'2020-02-26',1,[1]}...
    {'2020-02-28',2,[2]}...
    {'2020-03-02',1,[1]}...
    {'2020-03-05',1,[1]}...
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
%have commented out SVD fucntion which needs debugging 
[eyeData] = alignFace(expInfo, eyeData, behavioralData);

%%
%indexes trials by early and late 1st move and makes vectors for all Facemap ROIs
%plot graphs of pupil,whisking, paw and neural activity divided by early vs
%late trials
[earlyTrialsWhisk,lateTrialsWhisk,earlyTrialsPupil,lateTrialsPupil] = earlyVsLate(expInfo,behavioralData,neuralData,eyeData);

%% Pre-stimulus whisking analysis: compute linear fit and plot whisking +fit
%plot quantification of all early vs late trial pre-stim mean whisk & whisk slope
%needs to first run earlyVsLate function
prestimWhiskAnalysis(eyeData,earlyTrialsWhisk,lateTrialsWhisk,earlyTrialsPupil,lateTrialsPupil)
 
%% Rasters for all trials continously sorted (individual cells/ROIs) 
% use with 1 session loaded only 

%placing all trials in whichTrials  
%whichTrials=[];
%for t = 1:length(eyeData.eta.alignedFace{1}(:,1,1));
 %   whichTrials(end+1)=t;
%end

%selcting only trials which does not have pre-stim wheel movement.
[~, whichTrials] = selectCondition(expInfo, getUniqueContrasts(expInfo), behavioralData,...
    initTrialConditions('preStimMovement','quiescent','movementTime','late','specificRTs',[.1 3]));

%choosing cells correlating to whisking to plot only those
interp_whiskTraceNan = interp1(eyeData.timeAligned, eyeData.proc.face{1, 2}.motion, neuralData(1).respTimes); 
%remove the NaNs, but this is different for each session
interp_whiskTrace = interp_whiskTraceNan(~isnan(interp_whiskTraceNan));
nan_len = length(interp_whiskTraceNan)-length(interp_whiskTrace);
interp_whiskTrace = transpose(interp_whiskTrace);
interp_neural = neuralData.cellResps(nan_len+1:end,:);
whiskCells=[];
%start correlating from the time wihtout NaNss (approx sampling 50) in the interpolated whisktrace
for c=1: length(interp_neural(1,:));
    [r, p] = corrcoef(interp_whiskTrace, interp_neural(:,c));
    if p(1,2) < 0.05;
        whiskCells(end+1)= c;
    end
end
%choose your Raster plotting parameters:
%whichTrials=
whichCells = whiskCells;
k=1;
%plot raster of neural activity sorted by pre-stim whisking (200-0 ms)
figure;
rasterBrowser_whisking(expInfo, behavioralData, neuralData, whichCells, whichTrials, eyeData, 1)
figure;
rasterBrowser_facemap(expInfo, behavioralData, eyeData, 1:5, whichTrials, k)

%% Rasters for grouping trials (individual cells/ROIs)

%grouping trials by high vs low whisking (whisk-sorted) 
et = behavioralData.eventTimes;
wm = behavioralData.wheelMoves;
[relativeTimes,sortIdxWhisk] = sortTrialByWhiskgroup(whichTrials,eyeData,et,wm);
transpose(sortIdxWhisk);
%high shisking
whiskTrials{1} = (sortIdxWhisk(end-99:end));
%low
whiskTrials{2} = (sortIdxWhisk(1:100,1));

%filter cells that have significant activity during whiskTrials but not during nowhiskTrials 
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);
whiskCells=[];
for icell = 1:length(neuralData.eta.alignedResps{1,3})
    [p,h]=signrank(baselineResps(whiskTrials,icell),baselineResps(nowhiskTrials,icell));
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
plotFacemapRasters(expInfo, behavioralData, eyeData, 1:4, whichTrials,'byWhisk',k)

%% session raster for population neural activity compared to whisk-trace
%can change btween full session and specific smaller time window in the
%code 
plotSessionRaster(neuralData, eyeData)


