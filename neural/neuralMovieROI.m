for a = 1:size(videoSVD.motSVD_1,2)
    videoSVD.motSVD_int(:,a) = interp1(eyeData.timeAligned, videoSVD.motSVD_1(:,a)', neuralData.respTimes,'nearest');
end

%%

midpoint = floor(size(videoSVD.motSVD_int,1)/2);
Y = nanmean(neuralData.cellResps(:,isort1(1500:1700)),2);
X = videoSVD.motSVD_int;

X_test = X(1:midpoint,:);
Y_test = Y(1:midpoint,:);
X_train = X(midpoint:end,:);
Y_train = Y(midpoint:end,:);

%%

[B,dev,stats] = glmfit(X,Y);
%%
lambdas = [0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
[B,FitInfo] = lassoglm(X,Y,'normal','Lambda',lambdas);

%%

k = [0 1 10 100 1000 10000];
B = ridge(Y,X,k,1);

%%

for p = 1:500
    wMask(:,:,p) = videoSVD.motMask_reshape_1(:,:,p).*B(p,6);
end

figure;
% imfilt = imgaussfilt(sum(wMask,3),.5);
imagesc((sum(wMask,3)))
imagesc(imfilt)
colormap(BlueWhiteRed(100,.8))
caxis([-.001 .001])

%%

for a = 1:eyeData.timeAligned
    whichFrames = interp1(eyeData.timeAligned, 1:length(eyeData.timeAligned), neuralData.respTimes,'nearest');
end
%% construct movie matrix A (reshaped raw movie)
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
Abar = mean(A,1);
Asub = double(A) - Abar;

%% construct movie matrix Adiff (reshaped motion energy)
Adiff = zeros(length(whichFrames),160*214,'uint8');
for f = 1:length(whichFrames)
    try
    tmp1 = read(eyeData.veye,whichFrames(f));
    tmp2 = read(eyeData.veye,whichFrames(f)+1);
    ds1 = downsample(downsample(tmp1,3)',3)';
    ds2 = downsample(downsample(tmp2,3)',3)';
    dsdiff = abs(ds2 - ds1);
    Adiff(f,:) = reshape(dsdiff, [size(dsdiff,2)*size(dsdiff,1) 1]);
    catch
        Adiff(f,:) = nan(160*214,1);
    end
    prog = floor(100* (f/length(whichFrames)));
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end
% Anorm = double(A) / double(max(max(A)));
Adiff_bar = mean(Adiff,1);
Adiff_sub = double(Adiff) - Adiff_bar;


%% single cluster of cells (manual)
figure;
cellRange = 2700:2957;
Yz = zscore(nanmean(neuralData.cellResps(:,isort1(cellRange)),2));
pixWeights = (Yz'*Asub)/size(neuralData.cellResps,1);
pixMap = reshape(pixWeights,[size(ds,1) size(ds,2)]);

imagesc(imgaussfilt(pixMap,1));
colormap(BlueWhiteRed);
caxis([-4 4]);
axis off

%%
ROI_morph = zeros(size(ds,1), size(ds,2),size(neuralData.cellResps,2));
hop = 500;
for n = 1:size(neuralData.cellResps,2)
    try
        Yz = zscore(nanmean(neuralData.cellResps(:,isort1(n:n+hop)),2));
    catch
        Yz = zscore(nanmean(neuralData.cellResps(:,isort1(n:end)),2));
    end
pixWeights = (Yz'*Asub)/size(neuralData.cellResps,1);
pixMap = reshape(pixWeights,[size(ds,1) size(ds,2)]);
ROI_morph(:,:,n) = pixMap;
prog = floor(100* (n/size(neuralData.cellResps,2)));
    fprintf(1,'\b\b\b\b%3.0f%%',prog);
end

%%
fig1 = figure;
ax = subplot(1,4,1);
hl = line;
hr = line;
ht = line;
hb = line;
meanTrialActivity = squeeze(nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,sortIdx),1))';
meanPopActivity = squeeze(nanmean(neuralData.eta.alignedResps{1}(plotTrials,:,sortIdx),3))';
imagesc(eventWindow,1:length(normPSTH),meanTrialActivity);
colormap(ax,flipud(gray));
hold on;
line([0 0],[1 size(ROI_morph,3)],'LineStyle','--','Color','k')
caxis([0.05 .5])
xlim([-.5 2]);
set(gca,'tickdir','out')
box off

ax1 = subplot(1,4,[2 4]);
colormap(ax1,BlueWhiteRed);
caxis([-5 5]);
t = 1;

max_t = size(ROI_morph,3)-hop;
while t <= max_t
    
    subplot(ax);
    delete(hl);
    delete(hr);
    delete(ht);
    delete(hb);
    hl = line([-.5 -.5],[t t+hop],'LineStyle','-','Color','r','LineWidth',2);
    hr = line([2 2],[t t+hop],'LineStyle','-','Color','r','LineWidth',2);
    ht = line([-.5 2],[t t],'LineStyle','-','Color','r','LineWidth',2);
    hb = line([-.5 2],[t+hop t+hop],'LineStyle','-','Color','r','LineWidth',2);

    subplot(ax1);
    cla;
    imagesc(imgaussfilt(ROI_morph(:,:,t),1));
    hold on;
    axis off
    caxis([-5 5]);
    
    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig1, 'CurrentKey'), 'uparrow')
        t = max(1, t - 5);
    elseif was_a_key && strcmp(get(fig1, 'CurrentKey'), 'downarrow')
        t = min(max_t, t + 5);
    elseif was_a_key && strcmp(get(fig1, 'CurrentKey'), 'escape')
        close(fig1);
        break
    end
end

%%
clim = 150;
stds = [-3 -2 -1 0 1 2 3];
meanMap = reshape(Abar,[size(ds,1) size(ds,2)]);
figure;
imagesc(imgaussfilt(pixMap,1));
colormap(BlueWhiteRed);
caxis([-4 4]);
axis off
axis equal
% title(strcat({'face ROI (cell nos. '},num2str(cellRange(1)),'–',num2str(cellRange(end)),')'))
title('face ROI')


figure;
set(gcf,'position',[860 1280 2000 208])
for s = 1:length(stds)
    subplot(1,length(stds),s);
    
imagesc(meanMap + 10*stds(s)*pixMap)
colormap((gray));
caxis([-clim clim]);
axis off
axis tight
axis equal
title(strcat({'z = '},num2str(stds(s))))
end