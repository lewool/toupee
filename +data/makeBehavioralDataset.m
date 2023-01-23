% for m = 1:length(mouseList)
%     
%     mouseName = char(mouseList{m});
%     expDate = char(expList{m}{1});
%     expNum = expList{m}{2};
%     [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
%     fprintf('Finding %s...\n',expRef)
% 
%     sessDir = fullfile('G:\Analysis\Behavior',mouseName,expDate,num2str(expNum));
%     try 
%         cd(sessDir)
%         fprintf('Workspace exists!\n')
%     catch
%         mkdir(sessDir)
%         expInfo = initExpInfo({mouseList{m}},{expList{m}},hemList(m));
%         [expInfo, behavioralData] = processBehavior(expInfo);
%         cd(sessDir)
%         save('dataset.mat','behavioralData','expInfo')
%         fprintf('Workspace saved!\n')
%     end
%     
%     clearvars -except mouseList expList hemList m
% end
    
for m = 22:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Finding %s...\n',expRef)

    sessDir = fullfile('G:\Analysis\Behavior',mouseName,expDate,num2str(expNum));
    try 
        cd(sessDir)
        fprintf('Workspace exists - rerunning!\n')
        expInfo = initExpInfo({mouseList{m}},{expList{m}},hemList(m));
        [expInfo, behavioralData] = processBehavior(expInfo);
        cd(sessDir)
        save('dataset.mat','behavioralData','expInfo')
        fprintf('Workspace saved!\n')
    catch
        mkdir(sessDir)
        expInfo = initExpInfo({mouseList{m}},{expList{m}},hemList(m));
        [expInfo, behavioralData] = processBehavior(expInfo);
        cd(sessDir)
        save('dataset.mat','behavioralData','expInfo')
        fprintf('Workspace saved!\n')
    end
    
    clearvars -except mouseList expList hemList m
end