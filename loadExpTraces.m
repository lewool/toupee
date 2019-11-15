function [allFcell,ops] = loadExpTraces(mouseName, expDate, expSeries)
%Load the deconvolved spike data for curated ROIs across
%a series of experiments processed together in Suite2P

%% 
seriesFolder = strrep(num2str(expSeries),'  ',' ');
seriesFolder = strrep(num2str(seriesFolder),'  ',' ');
seriesFolder = strrep(num2str(seriesFolder),' ','_');

paths = data.dataPaths();
procDir = paths.local{1}{1};
file1Name = strcat('F_',mouseName,'_',expDate,'_plane',num2str(1),'_proc.mat');

load(fullfile(procDir,mouseName,expDate,seriesFolder,file1Name))
ops.numPlanes = dat.ops.nplanes;
ops.numChannels = dat.ops.nchannels;

allFcell = struct('spikes',[]);
for iPlane = 1:ops.numPlanes
    fileName = strcat('F_',mouseName,'_',expDate,'_plane',num2str(iPlane),'_proc.mat');
    load(fullfile(procDir,mouseName,expDate,seriesFolder,fileName))
    for iExp = 1:numel(dat.Fcell)

        spikesAll = dat.sp{1,iExp};

        %extract the dat.stat.iscell values as logicals to create an indexing vector
        iscellIdx = false(1,length(dat.stat));
        for iCell = 1:length(dat.stat)
            iscellIdx(iCell) = logical(dat.stat(iCell).iscell);
        end
        
        % extract the deconvolved spikes
        spikes = spikesAll(iscellIdx,:);
    
        planeTraces.spikes{1, iExp} = spikes;

        clearvars FcellAll FcellNeuAll Fcell FcellNeu neuropilCoeff FcellCorrected FcellCorrectedDFF FcellNeuDFF FcellNeuAllDFF spikes
    end
    %a single (corrected) cell trace is found in allPlaneTraces(iPlane).FcellCorrected{1,iExp}(iCell,:)
    allFcell(iPlane) = planeTraces;
    clear planeTraces;
end

