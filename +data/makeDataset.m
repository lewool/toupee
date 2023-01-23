for m = 1:length(mouseList)

    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Finding %s...\n',expRef)

    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    try 
        cd(sessDir)
        fprintf('Workspace exists!\n')
    catch
        try
            dataDir = fullfile('G:\Data\F',mouseName,expDate,num2str(expNum));
            cd(dataDir)
            mkdir(sessDir)
            expInfo = initExpInfo(mouseList(m),expList(m));
            [expInfo, neuralData, behavioralData] = processExperiment(expInfo);
            cd(sessDir)
            save('dataset.mat','behavioralData','expInfo','neuralData','-v7.3')
            fprintf('Workspace saved!\n')
        catch
            fprintf('No raw data exists!\n')
        end
    end
end
