clear all;
expInfo = initExpInfo({{'LEW031'}},{{'2020-02-03',1,[1]}});
expInfo.hemisphere = 1;
[expInfo, neuralData, behavioralData] = processExperiment(expInfo);
%%



%%
whichSessions = 1;
CRF = getEpochCRFs(expInfo(whichSessions), behavioralData(whichSessions),neuralData(whichSessions));

%%
maxY = [];
minY = [];
h5.fig = figure;
h5.fig.Position = [1250 160 1350 1170];
h5.ax = gobjects(4,4);

for ii = 1:16
    h5.ax(ii) = subplot(4,4,ii);
end

plotEpochCRFs(expInfo(whichSessions), neuralData(whichSessions), CRF, 'contraStim');
h1.fig = gcf;
for i = 1:4
    h1.ax(i) = subplot(1,4,i);
    xlabel('');
    if i == 1
        ylabel('contraStim')
    end
    maxY(end+1) = max(h1.ax(i).YAxis.Limits);
    minY(end+1) = min(h1.ax(i).YAxis.Limits);
end

plotEpochCRFs(expInfo(whichSessions), neuralData(whichSessions), CRF, 'ipsiStim');
h2.fig = gcf;
for i = 1:4
    h2.ax(i) = subplot(1,4,i);
    xlabel('');
    title('')
    if i ==1
        ylabel('ipsiStim')
    end
    maxY(end+1) = max(h2.ax(i).YAxis.Limits);
    minY(end+1) = min(h2.ax(i).YAxis.Limits);
end

plotEpochCRFs(expInfo(whichSessions), neuralData(whichSessions), CRF, 'contraMov');
h3.fig = gcf;
for i = 1:4
    h3.ax(i) = subplot(1,4,i);
    xlabel('');
    title('');
    if i == 1
        ylabel('contraMov')
    end
    maxY(end+1) = max(h3.ax(i).YAxis.Limits);
    minY(end+1) = min(h3.ax(i).YAxis.Limits);
end

plotEpochCRFs(expInfo(whichSessions), neuralData(whichSessions), CRF, 'ipsiMov');
h4.fig = gcf;
for i = 1:4
    h4.ax(i) = subplot(1,4,i);
    xlabel('');
    title('');
    if i ==1
        ylabel('ipsiMov')
    end
    maxY(end+1) = max(h4.ax(i).YAxis.Limits);
    minY(end+1) = min(h4.ax(i).YAxis.Limits);
end

h5.ax2 = gobjects(size(h5.ax));
for i = 1:4
    h5.ax2(i,1) = copyobj(h1.ax(i), h5.fig);
end

for i = 1:4
    h5.ax2(i,2) = copyobj(h2.ax(i), h5.fig);
end

for i = 1:4
    h5.ax2(i,3) = copyobj(h3.ax(i), h5.fig);
end

for i = 1:4
    h5.ax2(i,4) = copyobj(h4.ax(i), h5.fig);
end

for ii = 1:16
    h5.ax2(ii).Position = h5.ax(ii).Position;
    h5.ax2(ii).YAxis.Limits = [min(minY) max(maxY)];
end
delete(h5.ax);
close(h1.fig)
close(h2.fig)
close(h3.fig)
close(h4.fig)

%%
printfig(h5.fig,char(strcat(expInfo.mouseName,{' '},expInfo.expDate,{' CRFs'})))
close all;