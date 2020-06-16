function expFilePath = makeExpFilePath(expRef, expLog, dataFolder, fileType)
    fileType(1) = upper(fileType(1));
    suffix = strcat('_',fileType,'.mat');
    expFilePath = char(fullfile(dataFolder,expLog{1},expLog{2},expLog{3},strcat(expRef,suffix)));
end