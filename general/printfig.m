function printfig(gcf, name)

% get current working directory
cwd = cd;

% check wheter a filename has been provided; otherwise prompt
if nargin < 2
    prompt = {'Filename:'};
    dlgtitle = 'Print figure';
    dims = [1 35];
    definput = {''};
    figname = inputdlg(prompt,dlgtitle,dims,definput);
    name = figname{1};
end

% set the renderer
set(gcf,'renderer','Painters');
set(gcf, 'InvertHardCopy', 'off');
set(gcf, 'Color', 'w');

% go to the fig folder
cd('C:\Users\Wool\Desktop\tempFigs');

%print fig
print('-dpng', name, '-r300');

% return to the current working directory
cd(cwd);