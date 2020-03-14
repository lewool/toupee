function [allFcell,expInfo] = loadCellData(expInfo)
%Load the calcium for curated ROIs across
%a series of experiments processed together in Suite2P
% 24 Feb 2020 LEW added option to import processed data from python suite2p

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
file1Name = strcat('F_',mouseName,'_',expDate,'_plane',num2str(2),'_proc.mat');
file2Name = strcat('Fall.mat');

try
    load(fullfile(procDir,mouseName,expDate,seriesFolder,file1Name));
    expInfo.numPlanes = dat.ops.nplanes;
    expInfo.numChannels = dat.ops.nchannels;
    procType = 'matlab';
catch
    dat = load(fullfile(procDir,mouseName,expDate,seriesFolder,'suite2p','plane0',file2Name));
    expInfo.numPlanes = dat.ops.nplanes;
    expInfo.numChannels = dat.ops.nchannels;
    procType = 'python';
end

if strcmp(procType,'matlab') == 1
    for iPlane = 1:expInfo.numPlanes
        try
            fileName = strcat('F_',mouseName,'_',expDate,'_plane',num2str(iPlane),'_proc.mat');
            load(fullfile(procDir,mouseName,expDate,seriesFolder,fileName))
        catch
            fileName = strcat('F_',mouseName,'_',expDate,'_plane',num2str(iPlane),'.mat');
            warning(strcat('No processed file found for plane ',num2str(iPlane),'; loading default'));
            load(fullfile(procDir,mouseName,expDate,seriesFolder,fileName))
        end

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

elseif strcmp(procType,'python') == 1
    for iPlane = 1:expInfo.numPlanes
        folderName = strcat('plane',num2str(iPlane-1));
        dat = load(fullfile(procDir,mouseName,expDate,seriesFolder,'suite2p',folderName,file2Name));
        
        FcellAll = dat.F;
        FcellNeuAll = dat.Fneu;
        neuropilCoeff = dat.ops.neucoeff;
        spikesAll = dat.spks;
        
        %extract the dat.stat.iscell values as logicals to create an indexing vector
        iscellIdx = logical(dat.iscell(:,1));
        
        %reshape Fcell and FcellNeu to include only ROIs you selected manually during Suite2P
        Fcell = FcellAll(iscellIdx,:);
        FcellNeu = FcellNeuAll(iscellIdx,:);
        FcellCorrected = Fcell - (neuropilCoeff * FcellNeu);
        
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
        allFcell(iPlane) = planeTraces;
        clear planeTraces;
    end
end
