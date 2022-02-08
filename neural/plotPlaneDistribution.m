
for m = 1:length(mouseList)

mouseName = char(mouseList{m});
expDate = char(expList{m}{1});
expNum = expList{m}{2};
hemisphere = hemList(m);
[expRef, ~] = data.constructExpRef(mouseName,expDate,expNum);
fprintf('fitting %s...\n',expRef)
[behavioralData, expInfo, neuralData] = data.loadDataset(mouseName, expDate, expNum);
    
for p = 1:size(neuralData.allFcell,2)
    cellCounts(p) = size(neuralData.allFcell(p).Fcell,1);
    zlabels{p,1} = ones(cellCounts(p),1)*p;
end
all_Z = cell2mat(zlabels);
allStat = cat(1,neuralData.allFcell(:).stat);


if hemisphere < 0
    contraCells = find(neuralData.stats.bfcH(:,6));
    ipsiCells = find(neuralData.stats.bfcH(:,5));
    stimCCells = find(neuralData.stats.bfcH(:,3));
    stimICells = find(neuralData.stats.bfcH(:,2));
else
    contraCells = find(neuralData.stats.bfcH(:,5));
    ipsiCells = find(neuralData.stats.bfcH(:,6));
    stimCCells = find(neuralData.stats.bfcH(:,2));
    stimICells = find(neuralData.stats.bfcH(:,3));
end

[neuralData] = getSignificantActivity(expInfo, behavioralData, neuralData,0);

hitCells = find(neuralData.stats.bfcH(:,7));
missCells = find(neuralData.stats.bfcH(:,9));
dir_stimCells = find(neuralData.stats.bfcH(:,1) & (neuralData.stats.bfcH(:,6) | neuralData.stats.bfcH(:,5)));
dir_rewardCells = find(neuralData.stats.bfcH(:,7) & (neuralData.stats.bfcH(:,6) | neuralData.stats.bfcH(:,5)));
stim_rewardCells = find(neuralData.stats.bfcH(:,1) & neuralData.stats.bfcH(:,7));
tripleCells = find(neuralData.stats.bfcH(:,1) & neuralData.stats.bfcH(:,7) & (neuralData.stats.bfcH(:,6) | neuralData.stats.bfcH(:,5)));

%%
figure;
subplot(1,3,1)
for iPlane = 1:5
    z = iPlane * 40;
    im = (neuralData.allFcell(iPlane).meanImage);
    xImage = [0 512; 0 512];
    yImage = [0 0; 512 512]; 
    zImage = [z z; z z];
    surf(xImage,yImage,zImage,...    % Plot the surface
        'CData',im,...
        'FaceColor','texturemap');
    hold on
end
colormap(gray)
xlim([0 512]);
ylim([0 512]);
set(gca,'Zdir','reverse')
set(gca,'Ydir','reverse')
view(-30,20)
axis off

for iCell = 1:length(contraCells)
    xpix = allStat{contraCells(iCell), 1}.xpix;
    ypix = allStat{contraCells(iCell), 1}.ypix;
    zpix = ones(1,length(xpix))*all_Z(contraCells(iCell))*40;
    plot3(xpix,ypix,zpix,'.','Color',[1 0 0])
end
    
for iCell = 1:length(ipsiCells)
    xpix = allStat{ipsiCells(iCell), 1}.xpix;
    ypix = allStat{ipsiCells(iCell), 1}.ypix;
    zpix = ones(1,length(xpix))*all_Z(ipsiCells(iCell))*40;
    plot3(xpix,ypix,zpix,'.','Color',[0 .4 1])
end

subplot(1,3,2)
for iPlane = 1:5
    z = iPlane * 40;
    im = (neuralData.allFcell(iPlane).meanImage);
    xImage = [0 512; 0 512];
    yImage = [0 0; 512 512]; 
    zImage = [z z; z z];
    surf(xImage,yImage,zImage,...    % Plot the surface
        'CData',im,...
        'FaceColor','texturemap');
    hold on
end
colormap(gray)
xlim([0 512]);
ylim([0 512]);
set(gca,'Zdir','reverse')
set(gca,'Ydir','reverse')
view(-30,20)
axis off

for iCell = 1:length(stimCCells)
    xpix = allStat{stimCCells(iCell), 1}.xpix;
    ypix = allStat{stimCCells(iCell), 1}.ypix;
    zpix = ones(1,length(xpix))*all_Z(stimCCells(iCell))*40;
    plot3(xpix,ypix,zpix,'.','Color',[0 .4 1])
end

for iCell = 1:length(stimICells)
    xpix = allStat{stimICells(iCell), 1}.xpix;
    ypix = allStat{stimICells(iCell), 1}.ypix;
    zpix = ones(1,length(xpix))*all_Z(stimICells(iCell))*40;
    plot3(xpix,ypix,zpix,'.','Color',[1 0 0])
end

subplot(1,3,3)
for iPlane = 1:5
    z = iPlane * 40;
    im = (neuralData.allFcell(iPlane).meanImage);
    xImage = [0 512; 0 512];
    yImage = [0 0; 512 512]; 
    zImage = [z z; z z];
    surf(xImage,yImage,zImage,...    % Plot the surface
        'CData',im,...
        'FaceColor','texturemap');
    hold on
end
colormap(gray)
xlim([0 512]);
ylim([0 512]);
set(gca,'Zdir','reverse')
set(gca,'Ydir','reverse')
view(-30,20)
axis off

for iCell = 1:length(hitCells)
    xpix = allStat{hitCells(iCell), 1}.xpix;
    ypix = allStat{hitCells(iCell), 1}.ypix;
    zpix = ones(1,length(xpix))*all_Z(hitCells(iCell))*40;
    plot3(xpix,ypix,zpix,'.','Color',[0 .75 .1])
end

for iCell = 1:length(missCells)
    xpix = allStat{missCells(iCell), 1}.xpix;
    ypix = allStat{missCells(iCell), 1}.ypix;
    zpix = ones(1,length(xpix))*all_Z(missCells(iCell))*40;
    plot3(xpix,ypix,zpix,'.','Color',[.5 0 0])
end
%%
% 
% subplot(1,5,4)
% for iPlane = 1:5
%     z = iPlane * 40;
%     im = (neuralData.allFcell(iPlane).meanImage);
%     xImage = [0 512; 0 512];
%     yImage = [0 0; 512 512]; 
%     zImage = [z z; z z];
%     surf(xImage,yImage,zImage,...    % Plot the surface
%         'CData',im,...
%         'FaceColor','texturemap');
%     hold on
% end
% colormap(gray)
% xlim([0 512]);
% ylim([0 512]);
% set(gca,'Zdir','reverse')
% set(gca,'Ydir','reverse')
% view(-30,20)
% axis off
% 
% for iCell = 1:length(dir_stimCells)
%     xpix = allStat{dir_stimCells(iCell), 1}.xpix;
%     ypix = allStat{dir_stimCells(iCell), 1}.ypix;
%     zpix = ones(1,length(xpix))*all_Z(dir_stimCells(iCell))*40;
%     plot3(xpix,ypix,zpix,'.','Color',[1 .8 0])
% end
% 
% subplot(1,5,5)
% for iPlane = 1:5
%     z = iPlane * 40;
%     im = (neuralData.allFcell(iPlane).meanImage);
%     xImage = [0 512; 0 512];
%     yImage = [0 0; 512 512]; 
%     zImage = [z z; z z];
%     surf(xImage,yImage,zImage,...    % Plot the surface
%         'CData',im,...
%         'FaceColor','texturemap');
%     hold on
% end
% colormap(gray)
% xlim([0 512]);
% ylim([0 512]);
% set(gca,'Zdir','reverse')
% set(gca,'Ydir','reverse')
% view(-30,20)
% axis off
% 
% for iCell = 1:length(dir_rewardCells)
%     xpix = allStat{dir_rewardCells(iCell), 1}.xpix;
%     ypix = allStat{dir_rewardCells(iCell), 1}.ypix;
%     zpix = ones(1,length(xpix))*all_Z(dir_rewardCells(iCell))*40;
%     plot3(xpix,ypix,zpix,'.','Color',[1 .3 0])
% end


for z = 1:5
    total_dir_in_z = sum(all_Z(contraCells)==z) + sum(all_Z(ipsiCells)==z);
    contra_in_z(m,z) = sum(all_Z(contraCells)==z)/total_dir_in_z;
    ipsi_in_z(m,z) = sum(all_Z(ipsiCells)==z)/total_dir_in_z;
    
    total_stim_in_z = sum(all_Z(stimCCells)==z) + sum(all_Z(stimICells)==z);
    stimC_in_z(m,z) = sum(all_Z(stimCCells)==z)/total_stim_in_z;
    stimI_in_z(m,z) = sum(all_Z(stimICells)==z)/total_stim_in_z;
    
    total_rew_in_z = sum(all_Z(hitCells)==z) + sum(all_Z(missCells)==z);
    hit_in_z(m,z) = sum(all_Z(hitCells)==z)/total_rew_in_z;
    miss_in_z(m,z) = sum(all_Z(missCells)==z)/total_rew_in_z;
end

clearvars -except mouseList expList hemList contra_in_z ipsi_in_z stimC_in_z stimI_in_z hit_in_z miss_in_z
end
% figure;
% subplot(1,3,1)
% plot(contra_in_z,'Color',[0 .4 1],'LineWidth',2);
% hold on;
% plot(ipsi_in_z,'Color',[1 0 0],'LineWidth',2);
% 
% subplot(1,3,2)
% plot(stim_in_z,'Color',[0 .75 .1],'LineWidth',2);
% 
% subplot(1,3,3)
% plot(reward_in_z,'Color',[.5 0 .65],'LineWidth',2);