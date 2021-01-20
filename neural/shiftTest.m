
[baselineResps, stimResps, pmovResps, movResps, rewResps] = getEpochResps(neuralData.eta);

%% set up variables

% T is the total number of trials in the session (t indexes a single trial)
T = length(expInfo.block.events.endTrialTimes);

% N is the boundary for declaring the start and end of the central segment.
% While this is theoretically computed below, practically we want this to
% be ~19
N = 20;

% D is the central segment length
D = T - (2*N);

% X is a vector containing the block membership of each trial (-1 for left,
% +1 for right)

X = expInfo.block.events.highRewardSideValues(1:T);

for iCell = 931:931%1:length(baselineResps)
% Y is a vector of neural activity, where each element is a single value
% from each trial (e.g., mean baseline activity)
Y = baselineResps(:,iCell)';
% Y = rand(1,T);

%% central segment
% "We therefore extract a central segment of X, X[N:T-N], where N = (T - D)/2" 
% In matlab notation, this is X(N:T-N)

X0 = X(N:T-N);
Y0 = Y(N:T-N);

%% association value 
% "For real-valued series we might for example use the Pearson or Spearman
% correlation of X_t with Y_t over a segment of length D."

corrmat_0 = unique(corrcoef(X0, Y0));
V0_XY = corrmat_0(corrmat_0~=1);

shiftIndex = [-N+1:N-1];

% Vs_XY = zeros(1,length(shiftIndex));
for s = 1:length(shiftIndex)
    Ys = Y(N+shiftIndex(s):T-N+shiftIndex(s));
    corrmat = unique(corrcoef(X0, Ys));
    Vs_XY(iCell,s) = corrmat(corrmat~=1);
end

%% test statistic
% Compute a test statistic, m, that (1) determines whether a single shifted associated
% value (Vs) is equal to/greater than the nonshifted value (V0) for each
% shift (logical), then (2) sums up the logicals produced for each shift
% comparison.

% Basically, how many shifted comparisons generate an
% associate measure at least as big as the unshifted data?
a = 0.05;
m = length(find(Vs_XY(iCell,:) >= V0_XY))-1;
p = a*(2*N +1);
h(iCell) = m <= p || 2*N-m <=p;

end

%% plot your time series for an exampel cell
whichCell = 931;

figure;
set(gcf,'position',[38 558 1249 420]);
subplot(2,3,[1 2]);
plotD = fill([N T-N T-N N],[-1.2 -1.2 1.2 1.2],'b', 'LineStyle','none');
alpha(0.1);
hold on;
plot(X,'k');
ylim([-1.2 1.2]);
yticklabels({'L','','','','R'})
set(gca,'tickdir','out')
ylabel('block')
title('Time series')
xlim([1 T]);
box off

subplot(2,3,[4 5]);
plotD = fill([N T-N T-N N],[0 0 max(Y)*1.1 max(Y)*1.1],'b', 'LineStyle','none');
alpha(0.1);
hold on;
plot(Y,'k.');
xlim([1 T]);
ylim([0 max(Y)*1.1]);
xlabel('trial number')
ylabel('neural activity')
set(gca,'tickdir','out')
box off

subplot(2,3,[3 6])
line([0 0],[min(Vs_XY(whichCell,:))*1.1 max(Vs_XY(whichCell,:))*1.1],'Color',[.5 .5 .5],'LineStyle',':')
hold on
plot(shiftIndex,Vs_XY(whichCell,:),'k')
ylim([min(Vs_XY(whichCell,:))*1.1 max(Vs_XY(whichCell,:))*1.1]);
title('Shift test')
xlabel('shift')
ylabel('Pearson corr. coeff.')
set(gca,'tickdir','out')

text(-N*.9, min(Vs_XY(whichCell,:)),sprintf('m = %d \nh(0.05) = %d',m,h(whichCell)))
box off

