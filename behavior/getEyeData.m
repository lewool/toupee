function eyeData = getEyeData(expInfo)

for ex = 1:length(expInfo)
    paths = data.dataPaths();
    server = paths.server{1};
    eyefile = strcat(expInfo(ex).expDate,'_',num2str(expInfo(ex).expNum),'_',expInfo(ex).mouseName,'_eye.mat');
    eyeprocfile = strcat(expInfo(ex).expDate,'_',num2str(expInfo(ex).expNum),'_',expInfo(ex).mouseName,'_eye_proc.mat');
    
    try
        load(char(fullfile(server,expInfo(ex).mouseName,expInfo(ex).expDate,num2str(expInfo(ex).expNum),eyefile)));
        eyeData.time = extractfield(eyeLog.TriggerData,'Time');
    catch 
        error('No eye file found')
    end
    
    try
        proc = load(char(fullfile(server,expInfo(ex).mouseName,expInfo(ex).expDate,num2str(expInfo(ex).expNum),eyeprocfile)));
%         eyeData.proc = proc;
    catch
        warning('No processed eye file found')
    end
    
    %gather ROIs
    totalROIs = length(proc.rois);
    for r = 1:totalROIs
        if strcmp(proc.rois{r}.rtype  ,'pupil')
            eyeData.proc.pupil = proc.pupil{1, 1};
        elseif strcmp(proc.rois{r}.rtype  ,'motion SVD')
            fn = strcat('face',num2str(r-1));
            f = struct(...
                    'motion', proc.(matlab.lang.makeValidName(strcat('motion_',num2str(r-1)))),...
                    'motionSVD', proc.(matlab.lang.makeValidName(strcat('motSVD_',num2str(r-1)))));
            eyeData.proc.face{r-1} = f;
        end
    end   
end