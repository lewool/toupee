function CRF = getEpochCRFs(expInfo, behavioralData, neuralData)

for ex = 1:length(expInfo)
    
    %% get trial types and epoch responses (all cells, all trials)

    trialTypes = getTrialTypes(expInfo(ex), behavioralData(ex), 'late');
    whichTrials = trialTypes.intVar.all.contrast_direction;
    [baselineResps, stimIesps, pmovResps, movResps, rewResps] = getEpochResps(neuralData(ex).eta);

    %% multi-epoch CRF for different subpops of cells

    subpops = {'contraStim' 'ipsiStim' 'contraMov' 'ipsiMov' 'all'};

    for iS = 1:length(subpops)
        if strcmp(subpops{iS},'contraStim')
            if expInfo(ex).hemisphere > 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,2) == 1);
            elseif expInfo(ex).hemisphere < 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,3) == 1);
            end
        elseif strcmp(subpops{iS},'ipsiStim')
            if expInfo(ex).hemisphere > 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,3) == 1);
            elseif expInfo(ex).hemisphere < 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,2) == 1);
            end
        elseif strcmp(subpops{iS},'contraMov')
            if expInfo(ex).hemisphere > 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,5) == 1);
            elseif expInfo(ex).hemisphere < 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,6) == 1);
            end
        elseif strcmp(subpops{iS},'ipsiMov')
            if expInfo(ex).hemisphere > 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,6) == 1);
            elseif expInfo(ex).hemisphere < 0
                whichNeurons = find(neuralData(ex).stats.bfcH(:,5) == 1);
            end
        elseif strcmp(subpops{iS},'all')
            whichNeurons = (1:length(neuralData(ex).stats.bfcH))';
        end

        %initialize arrays
        bCRF_mC = zeros(size(whichTrials,1), size(whichNeurons,1));
        bCRF_mI = zeros(size(whichTrials,1), size(whichNeurons,1));
        sCRF_mC = zeros(size(whichTrials,1), size(whichNeurons,1));
        sCRF_mI = zeros(size(whichTrials,1), size(whichNeurons,1));
        pCRF_mC = zeros(size(whichTrials,1), size(whichNeurons,1));
        pCRF_mI = zeros(size(whichTrials,1), size(whichNeurons,1));
        mCRF_mC = zeros(size(whichTrials,1), size(whichNeurons,1));
        mCRF_mI = zeros(size(whichTrials,1), size(whichNeurons,1));

        %find neural response for all cells in each epoch, per contrast (mC vs mI)
        for c = 1:size(whichTrials,1)
            if expInfo(ex).hemisphere > 0
            trials_mC = whichTrials{c,1};
            trials_mI = whichTrials{c,2};
            elseif expInfo(ex).hemisphere < 0
                cflip = size(whichTrials,1) - c + 1;
                trials_mC = whichTrials{cflip,2};
                trials_mI = whichTrials{cflip,1};
            end

            bCRF_mC(c,:) = nanmean(baselineResps(trials_mC,whichNeurons),1);
            bCRF_mI(c,:) = nanmean(baselineResps(trials_mI,whichNeurons),1);

            sCRF_mC(c,:) = nanmean(stimIesps(trials_mC,whichNeurons),1);
            sCRF_mI(c,:) = nanmean(stimIesps(trials_mI,whichNeurons),1);

            pCRF_mC(c,:) = nanmean(pmovResps(trials_mC,whichNeurons),1);
            pCRF_mI(c,:) = nanmean(pmovResps(trials_mI,whichNeurons),1);

            mCRF_mC(c,:) = nanmean(movResps(trials_mC,whichNeurons),1);
            mCRF_mI(c,:) = nanmean(movResps(trials_mI,whichNeurons),1);

        end  
            
        if strcmp(subpops{iS},'contraStim')
            contraStim.mC.bCRF(:,ex) = nanmean(bCRF_mC,2);
            contraStim.mI.bCRF(:,ex) = nanmean(bCRF_mI,2);
            contraStim.mC.sCRF(:,ex) = nanmean(sCRF_mC,2);
            contraStim.mI.sCRF(:,ex) = nanmean(sCRF_mI,2);
            contraStim.mC.pCRF(:,ex) = nanmean(pCRF_mC,2);
            contraStim.mI.pCRF(:,ex) = nanmean(pCRF_mI,2);
            contraStim.mC.mCRF(:,ex) = nanmean(mCRF_mC,2);
            contraStim.mI.mCRF(:,ex) = nanmean(mCRF_mI,2);
            
        elseif strcmp(subpops{iS},'ipsiStim')
            ipsiStim.mC.bCRF(:,ex) = nanmean(bCRF_mC,2);
            ipsiStim.mI.bCRF(:,ex) = nanmean(bCRF_mI,2);
            ipsiStim.mC.sCRF(:,ex) = nanmean(sCRF_mC,2);
            ipsiStim.mI.sCRF(:,ex) = nanmean(sCRF_mI,2);
            ipsiStim.mC.pCRF(:,ex) = nanmean(pCRF_mC,2);
            ipsiStim.mI.pCRF(:,ex) = nanmean(pCRF_mI,2);
            ipsiStim.mC.mCRF(:,ex) = nanmean(mCRF_mC,2);
            ipsiStim.mI.mCRF(:,ex) = nanmean(mCRF_mI,2);

        elseif strcmp(subpops{iS},'contraMov')
            contraMov.mC.bCRF(:,ex) = nanmean(bCRF_mC,2);
            contraMov.mI.bCRF(:,ex) = nanmean(bCRF_mI,2);
            contraMov.mC.sCRF(:,ex) = nanmean(sCRF_mC,2);
            contraMov.mI.sCRF(:,ex) = nanmean(sCRF_mI,2);
            contraMov.mC.pCRF(:,ex) = nanmean(pCRF_mC,2);
            contraMov.mI.pCRF(:,ex) = nanmean(pCRF_mI,2);
            contraMov.mC.mCRF(:,ex) = nanmean(mCRF_mC,2);
            contraMov.mI.mCRF(:,ex) = nanmean(mCRF_mI,2);

        elseif strcmp(subpops{iS},'ipsiMov')
            ipsiMov.mC.bCRF(:,ex) = nanmean(bCRF_mC,2);
            ipsiMov.mI.bCRF(:,ex) = nanmean(bCRF_mI,2);
            ipsiMov.mC.sCRF(:,ex) = nanmean(sCRF_mC,2);
            ipsiMov.mI.sCRF(:,ex) = nanmean(sCRF_mI,2);
            ipsiMov.mC.pCRF(:,ex) = nanmean(pCRF_mC,2);
            ipsiMov.mI.pCRF(:,ex) = nanmean(pCRF_mI,2);
            ipsiMov.mC.mCRF(:,ex) = nanmean(mCRF_mC,2);
            ipsiMov.mI.mCRF(:,ex) = nanmean(mCRF_mI,2);

        elseif strcmp(subpops{iS},'all')
            all.mC.bCRF(:,ex) = nanmean(bCRF_mC,2);
            all.mI.bCRF(:,ex) = nanmean(bCRF_mI,2);
            all.mC.sCRF(:,ex) = nanmean(sCRF_mC,2);
            all.mI.sCRF(:,ex) = nanmean(sCRF_mI,2);
            all.mC.pCRF(:,ex) = nanmean(pCRF_mC,2);
            all.mI.pCRF(:,ex) = nanmean(pCRF_mI,2);
            all.mC.mCRF(:,ex) = nanmean(mCRF_mC,2);
            all.mI.mCRF(:,ex) = nanmean(mCRF_mI,2);
        end
    end
end

CRF = struct('contraStim',contraStim,'ipsiStim',ipsiStim,'contraMov',contraMov,'ipsiMov',ipsiMov,'all',all);
