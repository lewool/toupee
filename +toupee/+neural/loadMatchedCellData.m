function [allFcell, expInfo] = loadMatchedCellData(expInfo)
%Load the calcium for curated ROIs across
%a series of experiments processed separately, then
%cell-matched using ROImatchPub

%% GO THRU SESSIONS AND EXTRACT ONLY MATCHED CELLS

paths = data.dataPaths();
procDir = paths.local{1}{1};

listing = dir(fullfile(procDir,expInfo(1).subject,'matchedCells'));
numPlanes = length(listing) - 2;

for m = 1:length(expInfo)

    for iPlane = 1:numPlanes

        %load roiMatchData struct for the plane
        folderName = strcat('plane',num2str(iPlane));
        load(fullfile(procDir,expInfo(m).subject,'matchedCells',folderName,'matchFile.mat'));

        %load processed data for the session/plane
        dat = load(roiMatchData.allRois{m});
        
        %load expInfo
        expInfo(m).numPlanes = dat.ops.nplanes;
        expInfo(m).numChannels = dat.ops.nchannels;
        
        %load activity from the session/plane
        FcellAll = dat.F;
        FcellNeuAll = dat.Fneu;
        neuropilCoeff = dat.ops.neucoeff;
        spikesAll = dat.spks;

        %extract the session's iscell values to create an indexing vector
        iscellIdx = logical(dat.iscell(:,1));

        %reshape Fcell and FcellNeu to include only ROIs you selected
        %manually in the session
        FcellInd = FcellAll(iscellIdx,:);
        FcellNeuInd = FcellNeuAll(iscellIdx,:);

        %now extract only the cells that passed through cell-matching
        Fcell = FcellInd(roiMatchData.allSessionMapping(:,m),:);
        FcellNeu = FcellNeuInd(roiMatchData.allSessionMapping(:,m),:);

        %Neuropil correction
        FcellCorrected = Fcell - (neuropilCoeff * FcellNeu);

        % compute dF/F for each ROI trace/neuropil
        for j = 1:size(Fcell,1)
            FcellCorrectedDFF(j,:) = dff(FcellCorrected(j,:),0.1,20,30,numPlanes);
            FcellNeuDFF(j,:) = dff(FcellNeu(j,:),0.1,20,30,numPlanes);
        end

        % dF/F for ALL neuropil ROIs (NOT curated, NOT matched)
        for k = 1:size(FcellNeuAll,1)
            FcellNeuAllDFF(k,:) = dff(FcellNeuAll(k,:),0.1,20,30,numPlanes);
        end

        % extract the deconvolved spikes as well
        spikesInd = spikesAll(iscellIdx,:);
        spikes = spikesInd(roiMatchData.allSessionMapping(:,m),:);

        planeTraces.Fcell = Fcell;
        planeTraces.FcellNeu = FcellNeu;
        planeTraces.neuropilCoefficient = neuropilCoeff;
        planeTraces.FcellCorrected = FcellCorrected;
        planeTraces.FcellCorrectedDFF = FcellCorrectedDFF;
        planeTraces.FcellNeuDFF = FcellNeuDFF;
        planeTraces.FcellNeuAvg = mean(FcellNeuAll);
        planeTraces.FcellNeuAvgDFF = mean(FcellNeuAllDFF);
        planeTraces.spikes = spikes;

        clearvars FcellAll FcellNeuAll Fcell FcellNeu neuropilCoeff FcellCorrected FcellCorrectedDFF FcellNeuDFF FcellNeuAllDFF spikes

        %a single (corrected) cell trace is found in allPlaneTraces(iPlane).FcellCorrected(iCell,:)
        allFcellInd(iPlane) = planeTraces;
        clear planeTraces;
    end
    allFcell{m} = allFcellInd;
end
