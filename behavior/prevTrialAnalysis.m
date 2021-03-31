
%get trials of all different conditions, and convert them from logicals to list of all trials that meet the condictions  
for iX = 1:length(expInfo) 
    %find lists of trials with diff conditions happening in previous trial 
    [prevEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',[.1 3],'pastMovementTime','early'));
    [prevEarly{iX}] = find(prevEarly{iX}==1);
    
    [prevEarlyEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',[.1 3],'pastMovementTime','early','movementTime','early'));
    [prevEarlyEarly{iX}] = find(prevEarlyEarly{iX}==1);
    
    [prevLate{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastMovementTime','late'));
    [prevLate{iX}] = find(prevLate{iX}==1);
    
     [prevLateEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastMovementTime','late','movementTime','early'));
    [prevLateEarly{iX}] = find(prevLateEarly{iX}==1);

    [prevMiss{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','incorrect'));
    [prevMiss{iX}] = find(prevMiss{iX}==1);
    
    [prevMissEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','incorrect','movementTime','early'));
    [prevMissEarly{iX}] = find(prevMissEarly{iX}==1);

    [prevHit{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','correct'));
    [prevHit{iX}] = find(prevHit{iX}==1);
    
    [prevHitEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','correct','movementTime','early'));
    [prevHitEarly{iX}] = find(prevHitEarly{iX}==1);
   
    [allTrials{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('movementTime','all','preStimMovement','quiescent',...
        'specificRTs',([.1 3])));
    [allTrials{iX}] = find(allTrials{iX}==1);
    
    [preStimActive{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','active',...
        'specificRTs',([.1 3]),'pastResponseType','correct'));
    [preStimActive{iX}] = find(preStimActive{iX}==1);
    
    [preStimActiveEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','active',...
        'specificRTs',([.1 3]),'pastResponseType','correct','movementTime','early'));
    [preStimActiveEarly{iX}] = find(preStimActiveEarly{iX}==1);
    
    [preStimQuies{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','correct'));
    [preStimQuies{iX}] = find(preStimQuies{iX}==1);
    
    [preStimQuiesEarly{iX}, ~] = selectCondition(expInfo(iX), [-1 -.5 -.12 -.05 0 .05 .12 .5 1],...
        behavioralData(iX), initTrialConditions('preStimMovement','quiescent',...
        'specificRTs',([.1 3]),'pastResponseType','correct','movementTime','early'));
    [preStimQuiesEarly{iX}] = find(preStimQuiesEarly{iX}==1);
    
end

%% Assign, test and plot prev Early vs prev Late  
prevTrial = struct('prevTrialCond',[],'condNo',[],'fraction',[]);
condList = {{'Early'} {'Late'}};
for iX = 1:length(expList)
    b=0;
    for a = 1:length(condList)  
        %calculate the percentage of early trials vs ALL TRIALS for conditions
        if a == 1
            prevTrial(b+iX).fraction = 100 .* length(prevEarlyEarly{iX}) ./ length(prevEarly{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==2 
            prevTrial(b+iX).fraction = 100 .*length(prevLateEarly{iX}) ./ length(prevLate{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        end
        b=b+length(expInfo);
    end
end


%% VERSION 1 interactions : between condiction percentage early trials
%set up structure called prevTrials that will be called when plotting
prevTrial = struct('prevTrialCond',[],'condNo',[],'fraction',[]);
condList = {{'Early-hit'} {'Late-hit'} {'Early-miss'} {'Late-miss'}};

for iX = 1:length(expList)
    b=0;
    for a = 1:length(condList)  
        %calculate the percentage of early trials vs ALL TRIALS for conditions
        if a == 1
            prevTrial(b+iX).fraction = 100 .* length(prevEarlyHit{iX}) ./ length(allTrials{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==2 
            prevTrial(b+iX).fraction = 100 .*length(prevLateHit{iX}) ./ length(allTrials{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==3 
            prevTrial(b+iX).fraction = 100 .*length(prevEarlyMiss{iX}) ./ length(allTrials{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==4 
            prevTrial(b+iX).fraction = 100 .*length(prevLateMiss{iX}) ./ length(allTrials{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;            
        end
        b=b+length(expInfo);
    end
end


%% plot the trabsposed & extracted vector fields in a beeswarm plot  
figure;

beeswarm(extractfield(prevTrial,'condNo')', extractfield(prevTrial,'fraction')',...
    'sort_style','hex','dot_size',1,'overlay_style','sd','colormap','winter');
hold on
xlim([0 5])
ylim([0 30])
ylabel('%  Early trials of all trials','Fontsize',14)
%add significance lines!

%% do parametric stats test, i.e  ANOVA / glm
[p,tbl,stats] = anova1(extractfield(prevTrial,'fraction'),extractfield(prevTrial,'condNo'));
[c,m,h,nms]=multcompare(stats);

glm=fitlm(extractfield(prevTrial,'fraction'),extractfield(prevTrial,'condNo'));
res=glm.Residuals.raw;
qqplot(res)
%NO normality of residuals!

%% try non-parametric test 
[p,tbl,stats] = kruskalwallis(extractfield(prevTrial,'fraction'),extractfield(prevTrial,'condNo'));
c = multcompare(stats);

%% VERSION 2 interactions:   Wïthin conditionpercentages of early trials 
%set up structure called prevTrials that will be called when plotting
prevTrial = struct('prevTrialCond',[],'condNo',[],'fraction',[]);
condList = {{'Early-hit'} {'Late-hit'} {'Early-miss'} {'Late-miss'}};

for iX = 1:length(expList)
    b=0;
    for a = 1:length(condList)  
        %calculate the percentage of early trials vs ALL TRIALS for conditions
        if a == 1
            prevTrial(b+iX).fraction = 100 .* length(prevEarlyHitEarly{iX}) ./ length(prevEarlyHit{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==2 
            prevTrial(b+iX).fraction = 100 .*length(prevLateHitEarly{iX}) ./ length(prevLateHit{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==3 
            prevTrial(b+iX).fraction = 100 .*length(prevEarlyMissEarly{iX}) ./ length(prevEarlyMiss{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;
        elseif a==4 
            prevTrial(b+iX).fraction = 100 .*length(prevLateMissEarly{iX}) ./ length(prevLateMiss{iX});
            prevTrial(b+iX).prevTrialCond = condList{a};
            prevTrial(b+iX).condNo = a;            
        end
        b=b+length(expInfo);
    end
end


%% plot the trabsposed & extracted vector fields in a beeswarm plot  
figure;

beeswarm(extractfield(prevTrial,'condNo')', extractfield(prevTrial,'fraction')',...
    'sort_style','hex','dot_size',1,'overlay_style','sd','colormap','winter');
hold on
xlim([0 5])
ylim([0 30])
ylabel('%  Early trials within condition','Fontsize',14)
%add significance lines!


%%  non-parametric test 
[p,tbl,stats] = kruskalwallis(extractfield(prevTrial,'fraction'),extractfield(prevTrial,'condNo'));
c = multcompare(stats);





