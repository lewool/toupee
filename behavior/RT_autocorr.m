
%set up your vector of RTs, e.g., 
RT_vector = behavioralData(1).wheelMoves.epochs(5).onsetTimes' - behavioralData(1).eventTimes(1).daqTime';

%determine how many trials back you want to look at, let's say 5: 
numlags = 5;
%do a cross-correlation (an autocorrelation actually) of RT_vector with itself, but shifted 1-5 trials: 
[corr_RT, lags] = xcorr(RT_vector, RT_vector,numlags, 'coeff');
%we only need lags into the past 
corr_RT(lags<=0) = [];
lags(lags<=0) = [];
plot(lags,corr_RT,'ko')
xlabel('Past trials')
ylabel('Correlation coefficient')
title('Influence of previous reaction time')
%%
%do same for whik-RT corr 
whisk_vector = (eyeData(1).eta.alignedFace(:,95:101,2));
 
[variablesCorr, lags] = xcorr(RT_vector,whisk_vector, numlags, 'coeff');
variablesCorr(lags<=0) = [];
lags(lags<=0) = [];

%%
%plot
figure;
plot(lags,corr_RT,'ko')
xlabel('Past trials')
ylabel('Correlation coefficient')
title('Influence of previous reaction time')
hold on;
%plot(lags,variablesCorr,'ko')
%%
figure;
for t=2:length(RT_vector(:,1,1));
    plot(RT_vector(t-1,:,:),RT_vector(t,:,:),'ko')
    hold on
    xlim(1,20)
    ylim(1,20)
end
 
%%
%you can replace RT_vector with a vector of 0s/1s representing early/late
%and it should also work?
%and, instead of an autocorrelation, you can do the same correlation but with two different vectors, for example stimulus contrast or prestimulus whisking, move direction, or reward outcome, etc
variable1 = mean(eyeData(1).eta.alignFace(:,95:101));
variable2 = behavioralData(1).wheelMoves.epochs(5).onsetTimes' - behavioralData(1).eventTimes(1).daqTime';
 
[variablesCorr, lags] = xcorr(variable1, variable2, numlags, 'coeff');
variablesCorr(lags<=0) = [];
lags(lags<=0) = [];
plot(lags,variablesCorr,'ko')