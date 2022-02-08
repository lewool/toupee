function plotWhiskingVsVelocity(expInfo, behavioralData, eyeData, whichTrials)

et = behavioralData.eventTimes;
wm = behavioralData.wheelMoves;

[~, ~, meanPreStimWhisk] = sortTrialsByWhisking(whichTrials, eyeData, et, wm);

figure;
set(gcf,'position',[1490 1320 740 320])
subplot(1,2,1)
scatter(meanPreStimWhisk,abs(wm.epochs(5).peakVel(whichTrials)),10,'k','filled')
title(strcat(expInfo.mouseName,{' '},expInfo.expDate,{' (abs)'}))

subplot(1,2,2)
scatter(meanPreStimWhisk,(wm.epochs(5).peakVel(whichTrials)),10,'k','filled')
set(gca,'tickdir','out')
xlabel('mean pretrial whisking energy');
ylabel('peak wheel velocity');
title(strcat(expInfo.mouseName,{' '},expInfo.expDate,{' (signed)'}))