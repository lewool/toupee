subDim = ceil(sqrt(length(mouseList)));
figure(100);
useME = 0;
prog = 0;
for m = 14:length(mouseList)
    
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    hemisphere = hemList(m);
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    fprintf('fetching %s...',expRef)
    [behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    eyeData = getEyeData(expInfo);
    [eyeData] = alignFace(expInfo, eyeData, behavioralData);
    contrasts = getUniqueContrasts(expInfo);
    [~, lbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('highRewardSide','left','repeatType','random'));
    [~, rbTrials] = selectCondition(expInfo, contrasts, behavioralData, ...
        initTrialConditions('highRewardSide','right','repeatType','random'));
    
    %designate a baseline window
    Fs = .02;
    stim_eventIdx = find(eyeData.eta.eventWindow == 0);
    stim_preTime = [-0.5 0] / Fs;
    baselineIdx = stim_eventIdx + stim_preTime(1) : stim_eventIdx - 1;
    

    whichFrames = baselineIdx;
    fprintf(1,'movie matrix: %3d%%...',prog);
    if useME ~= 1
        
        %%%%%%% GENERATE A FROM RAW PIXELS
        A = zeros(size(eyeData.eta.alignedFrames{1},1),160*214,'uint8');
        for t = 1:size(eyeData.eta.alignedFrames{1},1)
            for f = 1:length(whichFrames)
                try
                tmp = read(eyeData.veye,eyeData.eta.alignedFrames{1}(t, whichFrames(f)));
                ds = downsample(downsample(tmp,3)',3)';
                A(t,:) = mean(reshape(ds, [size(ds,2)*size(ds,1) 1]),2);
                catch
                    A(t,:) = nan(size(ds,2)*size(ds,1),1);
                end
            end
            prog = floor(100* (t/size(eyeData.eta.alignedFrames{1},1)));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        
    else
        
        %%%%%%% GENERATE A FROM MOTION ENERGY
        A = zeros(size(eyeData.eta.alignedFrames{1},1),160*214,'uint8');
        for t = 1:size(eyeData.eta.alignedFrames{1},1)
            for f = 1:length(whichFrames)
                try
                tmp1 = read(eyeData.veye,eyeData.eta.alignedFrames{1}(t, whichFrames(f)));
                tmp2 = read(eyeData.veye,eyeData.eta.alignedFrames{1}(t, whichFrames(f))-1);
                diffframes = abs(tmp2-tmp1);
                ds = downsample(downsample(diffframes,3)',3)';
                A(t,:) = mean(reshape(ds, [size(ds,2)*size(ds,1) 1]),2);
                catch
                    A(t,:) = nan(size(ds,2)*size(ds,1),1);
                end
            end
            prog = floor(100* (t/size(eyeData.eta.alignedFrames{1},1)));
            fprintf(1,'\b\b\b\b%3.0f%%',prog);
        end
        
    end
    whichResps = double(A);
    
    if hemisphere > 0

        bC_trials = lbTrials;
        bI_trials = rbTrials;

    else
        bC_trials = rbTrials;
        bI_trials = lbTrials;

    end

    bC_resps = nanmean(whichResps(bC_trials,:),1)';
    bI_resps = nanmean(whichResps(bI_trials,:),1)';
    blockResp = bC_resps - bI_resps;


    %% generate pseudosessions
    prog = 0;
    fprintf(1,'movie matrix: %3d%%...',prog);
    nt = length(behavioralData.eventTimes(1).daqTime);
    for p = 1:1000
        b=zeros(1,nt);
        switches = cumsum(125+randi(100,1,10));
        for s = 1:length(switches)
            if s == 1
                b(1:switches(s)-1) = -1;
            elseif mod(s,2) == 1
                b(switches(s-1):switches(s)-1) = -1;
            elseif mod(s,2) == 0
                b(switches(s-1):switches(s)-1) = 1;
            end
        end
        b = b(1:nt);
        flip = randsample([-1, 1],1,true);
%         if b(1) ~= expInfo.block.paramsValues(1).firstHighSide
%             flip = -1;
%         else
%             flip = 1;
%         end
        b = flip*b;
        bC_trials_pseudo{p} = find(b < 0);
        bI_trials_pseudo{p} = find(b > 0);
        prog = floor(100* (p/1000));
        fprintf(1,'\b\b\b\b%3.0f%%',prog);
    end
    
    prog = 0;
    fprintf(1,'significant pixels: %3d%%...',prog);
    for p = 1:1000
        bC_resps_pseudo(:,p) = nanmean(whichResps(bC_trials_pseudo{p},:),1)';
        bI_resps_pseudo(:,p) = nanmean(whichResps(bI_trials_pseudo{p},:),1)';
        blockResp_pseudo(:,p) = bC_resps_pseudo(:,p) - bI_resps_pseudo(:,p);
        prog = floor(100* (p/1000));
        fprintf(1,'\b\b\b\b%3.0f%%',prog)
    end
    
    blockSig = zeros(1,length(blockResp));
    for c = 1:length(blockResp)
        UB = prctile(blockResp_pseudo(c,:),97.5);
        LB = prctile(blockResp_pseudo(c,:),2.5);
        if blockResp(c) < LB ||  blockResp(c) > UB
            blockSig(c) = 1;
        end
        [~, idx] = min(abs(sort(blockResp_pseudo(c,:)) - blockResp(c)));
        if idx > 500
            pValue{m}(c) = 2*((1000-idx)/1000);
        else
            pValue{m}(c) = 2*idx/1000;
        end
        sigrank{m}(c) = idx/1000;
    end

    propBS(m) = sum(blockSig)/length(blockResp);

    %% plot 'block value' of every cell for real vs pseudo sessions
    fprintf('plotting...')
%     set(gcf,'position',[32 80 2560 1556]);
    figure(100);
    mouseName = char(mouseList{m});
    expDate = char(expList{m}{1});
    expNum = expList{m}{2};
    [expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
    subplot(subDim,subDim,m)
    colors = [0.1 .7 .1; 1 .6 0];

    tmap{m} = reshape(mean(A,1),[160 214]);
    fmap{m} = imgaussfilt(reshape(sigrank{m},[160 214]),1.2);
    red1 = 1*cat(3, ones(size(tmap{m}))*colors(1,1), ones(size(tmap{m}))*colors(1,2), ones(size(tmap{m}))*colors(1,3));
    red2 = 0.5*cat(3, ones(size(tmap{m}))*colors(1,1), ones(size(tmap{m}))*colors(1,2), ones(size(tmap{m}))*colors(1,3));
    red3 = 0.25*cat(3, ones(size(tmap{m}))*colors(1,1), ones(size(tmap{m}))*colors(1,2), ones(size(tmap{m}))*colors(1,3));

    blue1 = 1*cat(3, ones(size(tmap{m}))*colors(2,1), ones(size(tmap{m}))*colors(2,2), ones(size(tmap{m}))*colors(2,3));
    blue2 = 0.5*cat(3, ones(size(tmap{m}))*colors(2,1), ones(size(tmap{m}))*colors(2,2), ones(size(tmap{m}))*colors(2,3));
    blue3 = 0.25*cat(3, ones(size(tmap{m}))*colors(2,1), ones(size(tmap{m}))*colors(2,2), ones(size(tmap{m}))*colors(2,3));

    imshow(tmap{m}(:,:,:,1));
    caxis([0 2]);
    hold on;

    h3 = image(red3);
    h2 = image(red2);
    h1 = image(red1);

    h4 = image(blue3);
    h5 = image(blue2);
    h6 = image(blue1);

    hold off
    set(h3, 'AlphaData', (fmap{m}>0.925))
    set(h2, 'AlphaData', (fmap{m}>0.95))
    set(h1, 'AlphaData', (fmap{m}>0.975))

    set(h4, 'AlphaData', (1-fmap{m})>0.925)
    set(h5, 'AlphaData', (1-fmap{m})>0.95)
    set(h6, 'AlphaData', (1-fmap{m})>0.975)
    axis off
    title(strcat(expRef),'Interpreter','none')
%%
    clearvars -except mouseList expList hemList pValue propBS subDim sigrank prog useME tmap fmap
    fprintf('done \n');
end













