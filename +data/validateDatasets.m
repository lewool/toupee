function [validatedMouseList, validatedExpList, idx] = validateDatasets(mouseList, expList)

vm = 1;
idx = false(1,length(mouseList));
for m = 1:length(mouseList)
    mouseName = mouseList{m};
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    
    
    
    try
        cd(sessDir)
        fileList = dir;
        if sum(contains({fileList.name},'data'))
            validatedMouseList{vm} = mouseList{m};
            validatedExpList{vm} = expList{m};
            vm = vm + 1;
            idx(m) = true;
        end
    catch 
    end
end