cd('G:\Workspaces');
for m = 1:length(expList)

mouseName = mouseList{1};
expDate = expList{m}{1};
expNum = expList{m}{2};

[expRef, expLog] = data.constructExpRef(mouseName,expDate,expNum);

sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
if ~isdir(sessDir{1})
    mkdir(sessDir{1})
end
end

%%
cd('G:\Workspaces');
for m = 12:length(expList)
    tic;
    
    mouseName = mouseList{m};
    expDate = expList{m}{1};
    expNum = expList{m}{2};
    
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir{1})
    
%     if ~exist('dataset.mat')
        fprintf('Loading %s:\n', char(expRef))
        expInfo = initExpInfo(mouseList(m),expList(m));
        [expInfo, neuralData, behavioralData] = processExperiment(expInfo);

        fprintf('Saving %s...', char(expRef))
        sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
        cd(sessDir{1})
        save dataset.mat -v7.3
        fprintf('done\n')   
%     end
    clearvars -except mouseList expList
    toc;
end

% for m = 1:length(mouseList)
%     plotPsychometric({mouseList{m}},{expList{m}});
%     printfig(gcf,strcat(char(mouseList{m}),'_',char(expList{m}{1}),'_psycho'));
%     close all
% end