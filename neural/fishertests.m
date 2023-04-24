allStrengths = {};
for m = 1:length(mouseList)
    
    %load data
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('Loading %s...\n',expRef)  
    sessDir = fullfile('G:\Workspaces',mouseName,expDate,num2str(expNum));
    cd(sessDir)
%     load('dataset.mat')
    load('kernelAnalysis.mat')

    propFeatures(m,:) = sum(kernelAnalysis.cellFeatureStrength(kernelAnalysis.maxEV > .01,:) > .01)./(sum(kernelAnalysis.maxEV > .01))*100;
    allStrengths{m,1} = kernelAnalysis.cellFeatureStrength(kernelAnalysis.maxEV > .01,:) > .01;
end

as = cat(1,allStrengths{:});
as_cells = [1:length(as)]';
combos = nchoosek([1 2 3 4 5 6],2);
%% all cells together
for aa = 1:length(combos)
    a = sum(as(:,combos(aa,1)) == 1 & as(:,combos(aa,2)) == 1);
    b = sum(as(:,combos(aa,1)) ~= 1 & as(:,combos(aa,2)) == 1);
    c = sum(as(:,combos(aa,1)) == 1 & as(:,combos(aa,2)) ~= 1);
    d = sum(as(:,combos(aa,1)) ~= 1 & as(:,combos(aa,2)) ~= 1);
    x = [a, b; c, d];
    [h(aa,:), p(aa,:)] = fishertest(x);
end

%% sessions separately

for m = 1:length(allStrengths)
    for aa = 1:length(combos)
        a = sum(allStrengths{m}(:,combos(aa,1)) == 1 & allStrengths{m}(:,combos(aa,2)) == 1);
        b = sum(allStrengths{m}(:,combos(aa,1)) ~= 1 & allStrengths{m}(:,combos(aa,2)) == 1);
        c = sum(allStrengths{m}(:,combos(aa,1)) == 1 & allStrengths{m}(:,combos(aa,2)) ~= 1);
        d = sum(allStrengths{m}(:,combos(aa,1)) ~= 1 & allStrengths{m}(:,combos(aa,2)) ~= 1);
        x = [a, b; c, d];
        [h(aa,m), p(aa,m)] = fishertest(x);
    end
end

figure;
line([0 7],[.05 .05],'Color',[.5 .5 .5],'LineStyle','--')
line([0 7],[.01 .01],'Color',[.5 .5 .5],'LineStyle','--')
hold on;
bar(median(p(median(h,2) == 1,:),2))
errorbar(median(p(median(h,2) == 1,:),2),mad(p(median(h,2) == 1,:),1,2),'linestyle','none','Color','k')
% errorbar(mean(p(median(h,2) == 1,:),2),std(p(median(h,2) == 1,:),1,2)/sqrt(36),'linestyle','none','Color','k')

ylim([0 .1])
prettyPlot(gca)

%%
sets = [1 2 3 4 5 6];
a=1;
for s = 1:length(sets)
    tmp = nchoosek(sets,s);
    for t=1:size(tmp,1)
        all_combinations{a} = tmp(t,:);
        a = a+1;
    end
end
for a = 1:length(all_combinations)
    test = false(1,6);
    test(all_combinations{a}) = true;
    occupants{a} = find(sum(as == test,2) == 6);
end

