function cm = customColormap(minColor, maxColor, n)

if nargin<3; n = 101; end
    
R = linspace(minColor(1),maxColor(1),n);
G = linspace(minColor(2),maxColor(2),n);
B = linspace(minColor(3),maxColor(3),n);

cm = [R' G' B'];
    
colormap(cm);