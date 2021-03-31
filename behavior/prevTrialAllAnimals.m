%% Import the data
[~, ~, raw] = xlsread('C:\Users\Ella Svahn\Documents\eyedata\PrevTrialEffect_allAnimals.xlsx','Sheet1','A2:G123');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

% Create output variable
data = reshape([raw{:}],size(raw));
% Create table
PT = table;

% Allocate imported array to column variable names
PT.animalID = data(:,1);
PT.propEarly_followingLate = data(:,2);
PT.propEarly_followingEarly = data(:,3);
PT.propEarly_followingCorrect = data(:,4);
PT.propEarly_followingIncorrect = data(:,5);
PT.propEarly_followingActive = data(:,6);
PT.propEarly_followingQuiescent = data(:,7);
% Clear temporary variables
clearvars data raw R;
%%

figure;
hold on;
subplot(1,2,1);

%plot the prop(early) for each session, paired between categories
p = plot([1, 2],[PT.propEarly_followingCorrect,PT.propEarly_followingIncorrect],...
    '-', 'LineWidth', 2,'Color',[.5 .5 .5],'MarkerSize', 10);
%transparency
%p.Color(4) = .5;
%don't erase 
hold on

%set x axis limits
xlim([.75 2.25])
%turn the box off
box off
%set tick direction
set(gca,'tickdir','out', 'Fontsize', 12)
%label only 1 and 2 on the x axis...
xticks([1 2])
%...and rename them
set(gca, 'XTickLabels', {'Rewarded', 'Non-rewarded'})
%label axes
xlabel('Previous trial condition')
ylabel('P(impulsive moves)')

subplot(1,2,2);

p = plot([1, 2],[PT.propEarly_followingEarly,PT.propEarly_followingLate], '-', 'LineWidth', 2, 'Color',[.5 .5 .5],'MarkerSize', 10);
%p.Color(4) = .5;
hold on

xlim([.75 2.25])
box off
set(gca,'tickdir','out','Fontsize', 12)
xticks([1 2])
set(gca, 'XTickLabels', {'Early', 'Late'})
xlabel('Previous trial condition')


%%
figure;
p = plot([1, 2],[PT.propEarly_followingActive,PT.propEarly_followingQuiescent], '-', 'LineWidth', 2, 'Color',[.5 .5 .5],'MarkerSize', 10);
%p.Color(4) = .5;
hold on

xlim([.75 2.25])
box off
set(gca,'tickdir','out','Fontsize', 12)
xticks([1 2])
set(gca, 'XTickLabels', {'Active pre-stim', 'Quiescent pre-stim'})
ylabel('P(impulsive moves)')

%xlabel('Trial condition')


%% do a Two-sample Kolmogorov-Smirnov test
disp('movetime')
[h,p] = kstest2(PT.propEarly_followingEarly,PT.propEarly_followingLate)
disp('reward')
[h,p] = kstest2(PT.propEarly_followingCorrect,PT.propEarly_followingIncorrect)
disp('pre-stimmove')
[h,p] = kstest2(PT.propEarly_followingActive,PT.propEarly_followingQuiescent)


%% beeswarm plot 
figure;

%beeswarm(extractfield(prevTrial,'condNo')', extractfield(prevTrial,'fraction')',...
 %   'sort_style','hex','dot_size',1,'overlay_style','sd','colormap','winter');

beeswarm(PT-animalID,[PT.propEarly_followingEarly,PT.propEarly_followingLate]',...
   'sort_style','hex','dot_size',1,'overlay_style','sd','colormap','winter');
hold on
%xlim([0 5])
%ylim([0 30])
ylabel('Fraction of early trials of all trials','Fontsize',14)
