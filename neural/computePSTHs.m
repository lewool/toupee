function [meanPSTH, semPSTH, rasters] = computePSTHs(cellArray,trialsArray)

% Cell array is a 2D trials x window array of responses for a single cell

for i = 1:size(trialsArray,1)
    try
        for j = 1:size(trialsArray,2)
            ti = trialsArray{i,j};

            %compute mean+sem responses
            rasters{i,j} = cellArray(ti,:);
            meanPSTH{i,j} = nanmean(cellArray(ti,:),1);
            semPSTH{i,j} = nanstd(cellArray(ti,:))/sqrt(length(ti));
        end
    catch
        for j = 1
            ti = trialsArray;

            %compute mean+sem responses
            rasters{i,j} = cellArray(ti,:);
            meanPSTH{i,j} = nanmean(cellArray(ti,:),1);
            semPSTH{i,j} = nanstd(cellArray(ti,:))/sqrt(length(ti));
        end
    end
end


% %contrasts
% allColors = [0 0 .5; 0 0 1; 0 .4 1; .6 .8 1; .75 .75 .75; .8 .45 .45; 1 0 0; .5 0 0; .25 0 0];
% 
% zeroIdx = find(contrasts == 0);
% walkup = length(contrasts) - zeroIdx;
% walkback = zeroIdx - 1;
% zeroGray = find(allColors(:,1) == .75);
% curveColors = allColors(zeroGray-walkback:zeroGray + walkup,:);
% 
% %sides
% curveColors = [0 .4 1; .75 .75 .75; 1 0 0];
% 
% %outcomes
% curveColors = [.1 .7 .1; .75 0 0];

