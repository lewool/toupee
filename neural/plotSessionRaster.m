function plotSessionRaster(neuralData, eyeData)

%% fetch cell responses & normalize
cellResps = neuralData.cellResps(4:end-3,:)';
respTimes = neuralData.respTimes(4:end-3);

%normalize each cell to 0-1
cmin = min(cellResps,[],2);
cb = cellResps - cmin;
cbmax = max(cb,[],2);
normCellResps = cb ./ cbmax;

%% fetch whisk data for the session 
whiskTrace = interp1(...
    eyeData.timeAligned, ...
    eyeData.proc.face{1, 2}.motion, ...
    neuralData.respTimes);
whiskTrace = whiskTrace(4:end-3);

%% compute whisking modulation index for each cell
%buffer NaNs at the start of whisk trace; ignore these
startCorrAt = 1 + max(find(isnan(whiskTrace)));

%compute the correlation coeff for each cell to the whisk trace
for iCell = 1:size(cellResps,1)
    [R,p] = corrcoef(cellResps(iCell,startCorrAt:end),whiskTrace(startCorrAt:end));
    whiskMod_Rvalues(iCell) = double(R(1,2));
    whiskMod_pValues(iCell) = double(p(1,2));
end

<<<<<<< HEAD
[~, sortIdx] = sort(whiskMod_Rvalues,'ascend');
=======
[~, sortIdx] = sort(whiskMod_Rvalues,'descend');
>>>>>>> master

%% plot all
figure;
set(gcf,'position',[57 155 1257 823])
subplot(8,1,[1:7])
imagesc(respTimes,1:size(normCellResps(sortIdx),1),normCellResps(sortIdx,:))
xlim([min(respTimes) max(respTimes)])
set(gca,'xtick',[])
set(gca,'ytick',[])
box off
set(gca,'tickdir','out')
<<<<<<< HEAD
ylabel('Neurons sorted by whisk correlation')
=======
ylabel('Neurons')
>>>>>>> master

colormap(flipud(gray));
caxis([0 .2])

subplot(8,1,8)
plot(respTimes,whiskTrace,'LineWidth',1.5,'Color',[1 0 1])
xlim([min(respTimes) max(respTimes)])
box off
set(gca,'tickdir','out')
xlabel('Session time (s)')
set(gca,'ytick',[])
ylabel('Whisking')

%% plot zoom
<<<<<<< HEAD
range = [4000 4150];
=======
range = [500 620];
>>>>>>> master
figure;
set(gcf,'position',[57 155 1257 823])
subplot(8,1,[1:7])
imagesc(respTimes,1:size(normCellResps(sortIdx),1),normCellResps(sortIdx,:))
xlim(range)
set(gca,'xtick',[])
set(gca,'ytick',[])
box off
set(gca,'tickdir','out')
<<<<<<< HEAD
ylabel('Neurons sorted by whisk correlation')
=======
ylabel('Neurons')
>>>>>>> master
colormap(flipud(gray));
caxis([0 .2])

subplot(8,1,8)
plot(respTimes,whiskTrace,'LineWidth',1.5,'Color',[1 0 1])
xlim(range)
set(gca,'ytick',[])
box off
set(gca,'tickdir','out')
xlabel('Session time (s)')
ylabel('Whisking')


