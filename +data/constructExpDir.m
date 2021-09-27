function expDir = constructExpDir(expRef, expLog, dataFolder)
    expDir = char(fullfile(dataFolder,expLog{1},expLog{2},expLog{3}));
end