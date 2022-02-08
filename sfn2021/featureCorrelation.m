for m = 1:length(mouseList)
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...\n',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    [neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);
%     hitMissTest{m} = hitMissTest(behavioralData, neuralData, expInfo, hemisphere);
    [fit{m} correlation{m}] = featureScatterplots(expInfo, behavioralData,neuralData, hemList(m));
    
    clearvars -except mouseList expList hemList ...
    hitResps_contra hitResps_ipsi missResps_contra missResps_ipsi ...
    contraResps_hit ipsiResps_miss contraResps_miss ipsiResps_miss ...
    hitResps_contraS missResps_contraS hitResps_ipsiS missResps_ipsiS hitMissTest ...
    fit correlation

end

%%

for m = 1:24
    f_ds(m) = fit{m}.dir_side(1);
    f_dr(m) = fit{m}.dir_reward(1);
    f_sr(m) = fit{m}.side_reward(1);
    
    corr_ds(m) = correlation{m}.dir_side;
    corr_dr(m) = correlation{m}.dir_reward;
    corr_sr(m) = correlation{m}.side_reward;
end

%%
figure;
scatter(1+(rand(1,24)-.5)*.1.*ones(1,24), corr_ds);
hold on;
scatter(2+(rand(1,24)-.5)*.1.*ones(1,24), corr_dr);
scatter(3+(rand(1,24)-.5)*.1.*ones(1,24), corr_sr);