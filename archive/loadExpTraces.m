function [allFcell,expInfo] = loadExpTraces(expInfo)
%Load the calcium for curated ROIs across
%a series of experiments processed together in Suite2P

%% LOAD EXPINFO VARIABLES

mouseName = expInfo.mouseName;
expDate = expInfo.expDate;
expSeries = expInfo.expSeries;

%% 
seriesFolder = strrep(num2str(expSeries),'  ',' ');
seriesFolder = strrep(num2str(seriesFolder),'  ',' ');
seriesFolder = strrep(num2str(seriesFolder),' ','_');

paths = data.dataPaths();
procDir = paths.local{1}{1};
file1Name = strcat('F_',mouseName,'_',expDate,'_plane',num2str(1),'_proc.mat');

load(fullfile(procDir,mouseName,expDate,seriesFolder,file1Name))
expInfo.numPlanes = dat.ops.nplanes;
expInfo.numChannels = dat.ops.nchannels;

for iPlane = 1:expInfo.numPlanes
    fileName = strcat('F_',mouseName,'_',expDate,'_plane',num2str(iPlane),'_proc.mat');
    load(fullfile(procDir,mouseName,expDate,seriesFolder,fileName))
    for iExp = 1:numel(dat.Fcell)

        FcellAll = dat.Fcell{1,iExp};
        FcellNeuAll = dat.FcellNeu{1,iExp};
        neuropilCoeffAll = extractfield(dat.stat,'neuropilCoefficient')';
        spikesAll = dat.sp{1,iExp};

        %extract the dat.stat.iscell values as logicals to create an indexing vector
        iscellIdx = false(1,length(dat.stat));
        for iCell = 1:length(dat.stat)
            iscellIdx(iCell) = logical(dat.stat(iCell).iscell);
        end
        
         %reshape Fcell and FcellNeu to include only ROIs you selected manually during Suite2P
        Fcell = FcellAll(iscellIdx,:);
        FcellNeu = FcellNeuAll(iscellIdx,:);
        neuropilCoeff = neuropilCoeffAll(iscellIdx,:);

        % Compute corrected cell trace (after neuropil subtraction)
        % Fcorrected = Fcell - (FcellNeu*NeuCoeff)
        for i = 1:size(Fcell,1)
            FcellCorrected(i,:) = Fcell(i,:) - (neuropilCoeff(i) .* FcellNeu(i,:));
        end
        
        % compute dF/F for each ROI trace/neuropil, plus all neuropil
        for j = 1:size(Fcell,1)
            FcellCorrectedDFF(j,:) = dff(FcellCorrected(j,:),0.1,20,30,expInfo.numPlanes);
            FcellNeuDFF(j,:) = dff(FcellNeu(j,:),0.1,20,30,expInfo.numPlanes);
        end
        
        for k = 1:size(FcellNeuAll,1)
            FcellNeuAllDFF(k,:) = dff(FcellNeuAll(k,:),0.1,20,30,expInfo.numPlanes);
        end
        
        % extract the deconvolved spikes as well
        spikes = spikesAll(iscellIdx,:);
        
        planeTraces.Fcell{1, iExp} = Fcell;
        planeTraces.FcellNeu{1, iExp} = FcellNeu;
        planeTraces.neuropilCoefficient{1, iExp} = neuropilCoeff;
        planeTraces.FcellCorrected{1, iExp} = FcellCorrected;
        planeTraces.FcellCorrectedDFF{1, iExp} = FcellCorrectedDFF;
        planeTraces.FcellNeuDFF{1, iExp} = FcellNeuDFF;
        planeTraces.FcellNeuAvg{1, iExp} = mean(FcellNeuAll);
        planeTraces.FcellNeuAvgDFF{1, iExp} = mean(FcellNeuAllDFF);
        planeTraces.spikes{1, iExp} = spikes;

        clearvars FcellAll FcellNeuAll Fcell FcellNeu neuropilCoeff FcellCorrected FcellCorrectedDFF FcellNeuDFF FcellNeuAllDFF spikes
    end
    %a single (corrected) cell trace is found in allPlaneTraces(iPlane).FcellCorrected{1,iExp}(iCell,:)
    allFcell(iPlane) = planeTraces;
    clear planeTraces;
end

