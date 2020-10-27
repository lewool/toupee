function [spidx1, spidx2] = goToRasterSubplot(xDim, yDim, rasterDims, whichX, whichRaster, buffer, startYAfter)
    spidx1 = sub2ind([xDim totalY], whichX, Ystart);
    spidx2 = sub2ind([xDim totalY], whichX, Yend);
    
    
Ystart = startYAfter+sum(rasterDims(1:whichRaster))-rasterDims(whichRaster)+buffer*whichRaster;
Yend = startYAfter+buffer*whichRaster+sum(rasterDims(1:whichRaster));