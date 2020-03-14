[plane1_pts,plane2_pts] = cpselect(mat2gray(plane1.img),mat2gray(plane2.img),'Wait',true);
tform = fitgeotrans(plane1_pts,plane2_pts,'affine');
% vfs_registered = imwarp(masked_vfs,tform,'OutputView',imref);
plane_registered = imwarp(plane1.img,tform);



%% single cell contour

%retrieve cell coords
figure;
hold on;
for c = 1:length(stat)
    if iscell(c,1) == 1
x = double(stat{c}.xpix)';
y = double(stat{c}.ypix)';

contIdx = boundary(x,y,.8);

cellContour = [x(contIdx) y(contIdx)];
plot(cellContour(:,1),cellContour(:,2),'b')
    end
end