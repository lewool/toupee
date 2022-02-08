function cm = colormapThruWhite(minColor, maxColor, n, gamma)

if nargin<1; n = 100; end
if nargin<2; gamma = 0.5; end

    midColor = [1 1 1];
    
    Rtop = linspace(midColor(1),maxColor(1),n/2)*n;
    Gtop = linspace(midColor(2),maxColor(2),n/2)*n;
    Btop = linspace(midColor(3),maxColor(3),n/2)*n;
    
    Rbottom = linspace(minColor(1), midColor(1),n/2)*n;
    Gbottom = linspace(minColor(2), midColor(2),n/2)*n;
    Bbottom = linspace(minColor(3), midColor(3),n/2)*n;
    
    rgb = [...
        Rbottom' Gbottom' Bbottom'; ...
        Rtop(2:end)' Gtop(2:end)' Btop(2:end)'; ...
        ];
    
    cm = (rgb/(n)).^gamma;
    
    %OLD
    % convert maxColor to HSV
%     maxHSV = rgb2hsv(maxColor);
%     cm = hsv2rgb(...
%         [...
%         ones(1,n)*maxHSV(1);
%         linspace(0,1,n);
%         ones(1,n)*maxHSV(3)...
%         ]');
   
colormap(cm);