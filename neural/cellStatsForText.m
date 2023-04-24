for m = 1:length(mouseList)

    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Finding %s...\n',expRef)

    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir)
    mm = load(fullfile(sessDir,'dataset.mat'));

    behavioralData = mm.behavioralData;
    expInfo = mm.expInfo;
    expInfo.hemisphere = hemList(m);
    neuralData = mm.neuralData;
    neuralData = getSingleCellStats(expInfo, behavioralData, neuralData);
       
    all_hTest{m,1} = neuralData.stats.hTest;
    clearvars -except mouseList expList hemList m all_hTest
end

%%

allCells = cat(1,all_hTest{:});

length(allCells(sum(allCells(:,1),2)>0,:))