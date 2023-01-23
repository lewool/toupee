eyeData = getEyeData(expInfo);
[eyeData] = alignFace(expInfo, eyeData, behavioralData);
%%

whichFrames = interp1(eyeData.timeAligned, 1:length(eyeData.timeAligned), neuralData.respTimes,'nearest');
whichFrames(isnan(whichFrames)) = 1;

%%
A = zeros(length(whichFrames),160*214,'uint8');
for f = 1:length(whichFrames)
    try
    tmp = read(eyeData.veye,whichFrames(f));
    ds = downsample(downsample(tmp,3)',3)';
    A(f,:) = reshape(ds, [160*214 1]);
    catch
        A(f,:) = nan(size(160*214,1));
    end
    prog = floor(100* (f/length(whichFrames)));
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end
% Anorm = double(A) / double(max(max(A)));
% Abar = mean(A,1);
% Asub = double(A) - Abar;

%%

A = zeros(length(whichFrames),120*160,'uint8');
for f = 1:length(whichFrames)
    try
    tmp = read(eyeData.veye,whichFrames(f));
    ds = downsample(downsample(tmp,4)',4)';
    A(f,:) = reshape(ds, [120*160 1]);
    catch
        A(f,:) = nan(size(120*160,1));
    end
    prog = floor(100* (f/length(whichFrames)));
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end

%%
Adiff = [zeros(1, 120*160); abs(diff(A,1))];

%%

motionEnergy = motSVD_1(:,1);
for ETA = 1:4
for t = 1:size(eyeData.eta.alignedFrames{1},1)
    try
        alignedME{ETA}(t,:) = motionEnergy(eyeData.eta.alignedFrames{ETA}(t,:));
    catch
        alignedME{ETA}(t,:) = nan;
    end
end
end

%%
colors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
% figure;
contrasts = getUniqueContrasts(expInfo);
trialTypes = getTrialTypes(expInfo, behavioralData, 'early');
subplot(1,2,2)
for c = 1:9
    if c < 5
        whichTrials = trialTypes.intVar.all.contrast_direction{c,1};
        plot(contrasts(c),-mean(behavioralData.wheelMoves.epochs(5).peakVel(whichTrials)),'o','Color',colors(c,:))
        hold on
    elseif c > 5
        whichTrials = trialTypes.intVar.all.contrast_direction{c,2};
        plot(contrasts(c),-mean(behavioralData.wheelMoves.epochs(5).peakVel(whichTrials)),'o','Color',colors(c,:))
        hold on
    else
        wA = trialTypes.intVar.all.contrast_direction{c,1};
        wB = trialTypes.intVar.all.contrast_direction{c,2};
        plot(contrasts(c),-mean(behavioralData.wheelMoves.epochs(5).peakVel(wA)),'o','Color',colors(c,:))
        plot(contrasts(c),-mean(behavioralData.wheelMoves.epochs(5).peakVel(wB)),'o','Color',colors(c,:))
    end

end
ylim([-120 120])


prettyPlot(gca)