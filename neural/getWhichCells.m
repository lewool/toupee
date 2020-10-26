function plotCells = getWhichCells(whichCells,neuralData)

% whichCells = 'leftStim'; %choose from 'pLabels' array
if strcmp(whichCells, 'all')
    plotCells = 1:size(neuralData.eta.alignedResps{1},3);
elseif strcmp(whichCells, 'advanceRight')
    plotCells = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'advanceMov')) > 0 & neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'rightMov')) > 0);
elseif strcmp(whichCells, 'advanceLeft')
    plotCells = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'advanceMov')) > 0 & neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,'leftMov')) > 0);
else
    plotCells = find(neuralData.stats.bfcH(:,strcmp(neuralData.stats.labels,whichCells)) > 0);
end