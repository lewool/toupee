percentSig = sum(neuralData.stats.bfcH,1)/length(neuralData.stats.bfcH);
lAD = sum(bfcH(:,9) .* bfcH(:,5))/length(neuralData.stats.bfcH);
rAD = sum(bfcH(:,9) .* bfcH(:,6))/length(neuralData.stats.bfcH);

figure;hold on
set(gcf,'position',[127 730 1320 248])
subplot(1,7,6)
hPieComponentHandles = pie([percentSig(1) 1-percentSig(1)]);
set(hPieComponentHandles(1),'FaceColor','k')
set(hPieComponentHandles(3),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')
xlabel('stim')

subplot(1,7,1)
hPieComponentHandles = pie([percentSig(4) 1-percentSig(4)]);
set(hPieComponentHandles(1),'FaceColor',[.7 0 .7])
set(hPieComponentHandles(3),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')

subplot(1,7,2)
hPieComponentHandles = pie([percentSig(5) percentSig(6) 1-percentSig(5)-percentSig(6)]);
set(hPieComponentHandles(1),'FaceColor',[0 .4 1])
set(hPieComponentHandles(3),'FaceColor',[1 0 0])
set(hPieComponentHandles(5),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')
set(hPieComponentHandles(5),'EdgeColor','none')

subplot(1,7,5)
hPieComponentHandles = pie([percentSig(7) 1-percentSig(7)]);
set(hPieComponentHandles(1),'FaceColor',[.2 .8 0])
set(hPieComponentHandles(3),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')


subplot(1,7,3)
hPieComponentHandles = pie([percentSig(9) 1-percentSig(9)]);
set(hPieComponentHandles(1),'FaceColor',[1 .8 0])
set(hPieComponentHandles(3),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')

subplot(1,7,4)
hPieComponentHandles = pie([lAD rAD 1-(lAD + rAD)]);
set(hPieComponentHandles(1),'FaceColor',[0 .4 1])
set(hPieComponentHandles(3),'FaceColor',[1 0 0])
set(hPieComponentHandles(5),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')
set(hPieComponentHandles(5),'EdgeColor','none')

subplot(1,7,7)
hPieComponentHandles = pie([percentSig(8) 1-percentSig(8)]);
set(hPieComponentHandles(1),'FaceColor',[1 .5 0])
set(hPieComponentHandles(3),'FaceColor',[.75 .75 .75])
set(hPieComponentHandles(1),'EdgeColor','none')
set(hPieComponentHandles(3),'EdgeColor','none')

