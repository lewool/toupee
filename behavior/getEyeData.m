function eyeData = getEyeData(expInfo)

for ex = 1:length(expInfo)
    server = expInfo(ex).server;
    eyefile = strcat(expInfo(ex).expDate,'_',num2str(expInfo(ex).expNum),'_',expInfo(ex).mouseName,'_eye.mat');
    eyeprocfile = strcat(expInfo(ex).expDate,'_',num2str(expInfo(ex).expNum),'_',expInfo(ex).mouseName,'_eye_proc.mat');
    
    try
        load(fullfile(server,expInfo(ex).mouseName,expInfo(ex).expDate,num2str(expInfo(ex).expNum),eyefile));
        eyeData.time = extractfield(eyeLog.TriggerData,'Time');
    catch 
        error('No eye file found')
    end
    
    try
        load(fullfile(server,expInfo(ex).mouseName,expInfo(ex).expDate,num2str(expInfo(ex).expNum),eyeprocfile));
        eyeData.proc = proc;
    catch
        warning('No processed eye file found')
    end
end