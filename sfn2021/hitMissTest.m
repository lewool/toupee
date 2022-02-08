function hitMissTest = hitMissTest(behavioralData, neuralData, expInfo, hemisphere)

ns = 1000;
a = 0.05;

for s = 1:ns
    trialTypes = getTrialTypes(expInfo, behavioralData, 'late');
    clc;
    fprintf('%d \n',s)
    [baselineResps, ~, ~, ~, ~, ~] = getEpochResps(neuralData.eta);

    whichETA = 1;

    hitCells = find(neuralData.stats.bfcH(:,7));
    missCells = find(neuralData.stats.bfcH(:,9));
    sL_trials = trialTypes.singleVar.side{1};
    sR_trials = trialTypes.singleVar.side{3};
    mL_trials = trialTypes.singleVar.direction{1};
    mR_trials = trialTypes.singleVar.direction{2};
    hit_trials = trialTypes.singleVar.outcome{1};
    miss_trials = trialTypes.singleVar.outcome{2};
    mLhit_trials = trialTypes.intVar.cb2D.outcome_direction{1,1};
    mRhit_trials = trialTypes.intVar.cb2D.outcome_direction{1,2};
    mLmiss_trials = trialTypes.intVar.cb2D.outcome_direction{2,1};
    mRmiss_trials = trialTypes.intVar.cb2D.outcome_direction{2,2};

    resps_mLhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLhit_trials,:,leftMovCells)))' - mean(baselineResps(mLhit_trials,leftMovCells),1)';
    resps_mRhitTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRhit_trials,:,leftMovCells)))' - mean(baselineResps(mRhit_trials,leftMovCells),1)';
    resps_mLhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLhit_trials,:,rightMovCells)))' - mean(baselineResps(mLhit_trials,rightMovCells),1)';
    resps_mRhitTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRhit_trials,:,rightMovCells)))' - mean(baselineResps(mRhit_trials,rightMovCells),1)';

    resps_mLmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLmiss_trials,:,leftMovCells)))' - mean(baselineResps(mLmiss_trials,leftMovCells),1)';
    resps_mRmissTrials_lCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRmiss_trials,:,leftMovCells)))' - mean(baselineResps(mRmiss_trials,leftMovCells),1)';
    resps_mLmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mLmiss_trials,:,rightMovCells)))' - mean(baselineResps(mLmiss_trials,rightMovCells),1)';
    resps_mRmissTrials_rCells = squeeze(mean(neuralData.eta.alignedResps{whichETA}(mRmiss_trials,:,rightMovCells)))' - mean(baselineResps(mRmiss_trials,rightMovCells),1)';

    if hemisphere > 0
        %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
        dHc = resps_mLhitTrials_lCells - resps_mRhitTrials_lCells;
        dHi = resps_mLhitTrials_rCells - resps_mRhitTrials_rCells;
        dMc = resps_mLmissTrials_lCells - resps_mRmissTrials_lCells;
        dMi = resps_mLmissTrials_rCells - resps_mRmissTrials_rCells;
    else
        %deltaResp = contraTrials - ipsiTrials (so contra responses are positive)
        dHc = resps_mRhitTrials_rCells - resps_mLhitTrials_rCells;
        dHi = resps_mRhitTrials_lCells - resps_mLhitTrials_lCells;
        dMc = resps_mRmissTrials_rCells - resps_mLmissTrials_rCells;
        dMi = resps_mRmissTrials_lCells - resps_mLmissTrials_lCells;
    end

    for t = 1:size(dHc,2)
        [~, h_hitMissTest_contra(s,t)] = signrank(dHc(:,t),dMc(:,t));
        [~, h_hitMissTest_ipsi(s,t)] = signrank(dHi(:,t),dMi(:,t));
    end
end

hitMissTest.hOverTimebin_contra = (1-sum(h_hitMissTest_contra)/ns) < a;
hitMissTest.hOverTimebin_ipsi = (1-sum(h_hitMissTest_ipsi)/ns) < a;
hitMissTest.alpha = a;
hitMissTest.iterations = ns;