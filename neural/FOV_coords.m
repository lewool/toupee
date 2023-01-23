function FOV_xyz = FOV_coords(expInfo, neuralData)
FOV_xyz = [];
spacebetweenplanes=60;

for iPlane = 1:size(neuralData.allFcell,2)

    cellCounts = numel(neuralData.allFcell(iPlane).stat);
    if iPlane == 1

        startCount = 1;
        endCount =  cellCounts;

    elseif iPlane > 1
        startCount = endCount +1;
        endCount =  endCount+cellCounts;  

    end

    for c =1:cellCounts
        FOV_xyz(c+startCount-1).x =mean(neuralData.allFcell(iPlane).stat{c,1}.xpix);
        FOV_xyz(c+startCount-1).y =mean(neuralData.allFcell(iPlane).stat{c,1}.ypix);
        FOV_xyz(c+startCount-1).z= spacebetweenplanes*iPlane;
    end
end
