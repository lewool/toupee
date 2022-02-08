function cm = RedWhiteGreen(n, gamma)

if nargin<1; n = 100; end
if nargin<2; gamma = 0.5; end

cm = ([0:n, ones(1,n)*n ; ...
      n*ones(1,n), n:-1:0; ...
      0:n, n-1:-1:0]' / n).^gamma;

  cm = flipud(cm);
colormap(cm);