function cm = rasterColor(maxColor, n)

if nargin<1; maxColor = [0 0 0]; end
if nargin<2; n = 101; end

if maxColor == [0 0 0]
    cm = flipud(gray(n));
else
    minColor = [1 1 1];
    
    R = linspace(minColor(1),maxColor(1),n);
    G = linspace(minColor(2),maxColor(2),n);
    B = linspace(minColor(3),maxColor(3),n);
    
    cm = [R' G' B'];
    
    %OLD
    % convert maxColor to HSV
%     maxHSV = rgb2hsv(maxColor);
%     cm = hsv2rgb(...
%         [...
%         ones(1,n)*maxHSV(1);
%         linspace(0,1,n);
%         ones(1,n)*maxHSV(3)...
%         ]');
    
end
colormap(cm);