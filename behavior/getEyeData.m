function eyeData = getEyeData(expInfo)

for ex = 1:length(expInfo)
    paths = data.dataPaths();
    server = paths.server{2};
    eyefile = strcat(expInfo(ex).expDate,'_',num2str(expInfo(ex).expNum),'_',expInfo(ex).mouseName,'_eye.mat');
    eyeprocfile = strcat(expInfo(ex).expDate,'_',num2str(expInfo(ex).expNum),'_',expInfo(ex).mouseName,'_eye_proc.mat');
    
    try
        warning off
        load(char(fullfile(server,expInfo(ex).mouseName,expInfo(ex).expDate,num2str(expInfo(ex).expNum),eyefile)));
        eyeData(ex).time = extractfield(eyeLog.TriggerData,'Time');
        warning on
    catch 
        warning on
        error('No eye file found')
    end
    
    try
        warning off
        proc = load(char(fullfile(server,expInfo(ex).mouseName,expInfo(ex).expDate,num2str(expInfo(ex).expNum),eyeprocfile)));
        warning on
%         eyeData.proc = proc;
    catch
        warning on
        warning('No processed eye file found')
    end
    
    if ~isequal(length(proc.pupil{1, 1}.area),length(eyeData(ex).time))
        eyeData(ex).time = eyeData(ex).time(1:length(proc.pupil{1, 1}.area));
    end
    
    %gather ROIs
    totalROIs = length(proc.rois);
    for r = 1:totalROIs
        if strcmp(proc.rois{r}.rtype  ,'pupil')
            eyeData(ex).proc.pupil = proc.pupil{1, 1};
            eyeData(ex).proc.pupil.area = zscore(eyeData(ex).proc.pupil.area);
        elseif strcmp(proc.rois{r}.rtype  ,'motion SVD')
            f = struct(...
                    'motion', proc.(matlab.lang.makeValidName(strcat('motion_',num2str(r-1)))),...
                    'motionSVD', proc.(matlab.lang.makeValidName(strcat('motSVD_',num2str(r-1)))));
                eyeData(ex).proc.face{r-1} = f;
                eyeData(ex).proc.face{r-1}.motion = zscore(eyeData(ex).proc.face{r-1}.motion);
        end
    end
    
    %align video timestams to Timeline using xcorr
    lickMotion = eyeData(ex).proc.face{2}.motion;
    lickport = diff(expInfo(ex).Timeline.rawDAQData(:,12));

    lickMotion_int = interp1(...
        eyeData(ex).time, ...
        lickMotion, ...
        expInfo(ex).Timeline.rawDAQTimestamps(1:end-1));
    
    lickMotion_int(isnan(lickMotion_int)) = 0;

    [c, lags] = xcorr(lickMotion_int, lickport);

    [~,i] = max(c);
    timeShift = lags(i)*0.001;

    eyeData(ex).timeAligned = eyeData(ex).time - timeShift;
    
    % plot to check alignment
    expRef = strcat({expInfo(ex).mouseName},{'_'},{expInfo(ex).expDate},{'_'},{num2str(expInfo(ex).expNum)});
    figure;
    set(gcf,'position',[44 635 1767 343]);
    p1 = plot(expInfo(ex).Timeline.rawDAQTimestamps(1:end-1),diff(expInfo(ex).Timeline.rawDAQData(:,12)),'Color',[.7 .7 .7]);
    p1.Color(4) = 1;
    ax = gca;
    title(strcat({'lick data alignment: '},{char(expRef)},{' (shift = '},{num2str(timeShift)},{' s)'}),'Interpreter','none')
    ax.TickDir = 'out';
    xlabel('time (seconds)')
    ylabel('optical beam events')
    xlim([0 100]);
    box off
    yyaxis right;
    p2 = plot(eyeData(ex).timeAligned, eyeData(ex).proc.face{2}.motion,'m');
    p2.Color(4) = .5;
    xlim([0 100]);
    ylabel('video ROI motion SVD')
    set(gca, 'YColor', 'm');
    box off
    
end