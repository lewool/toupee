

C = sum(importdata('2020-08-13_LEW040_window.tif'),3)/255;
% C = fliplr(C);
%%

FOV = importdata('2020-09-23_LEW040_FOV_1.7x.jpg');
figure;imagesc(FOV);
axis equal
axis off
%%

fig = figure; imagesc(C);
colormap((gray));
hold on
axis off
axis equal
t = 1; 
x = 914;
y = 590;
r = 290;
h0 = line; h1 = line; h2 = line;

while t < 10000000

    delete(h0);
    delete(h1);
    delete(h2);
   
    hold on;
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h0 = plot(xunit, yunit,'r');
    h1 = line([x x],[y-r y+r],'Color','r');
    h2 = line([x-r x+r],[y y],'Color','r');

    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
        x = x - 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
        x = x + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'a')
        r = r - 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'uparrow')
        y = y - 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'downarrow')
        y = y + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 's')
        r = r + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        coords = [x y r];
        break        
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        close(fig); break
    end
end
    
%%

fig = figure;
set(gcf,'position',[1256 440 1538 1186]);
subplot(1,3,1)
imagesc(fliplr(FOV));
axis equal
axis off
subplot(1,3,[2 3])
imagesc(C);
colormap((gray));
hold on
axis off
axis equal
% axis ij
t = 1; 
x = coords(1);
y = coords(2);
r = (coords(3)*0.4)/2;
h0 = line; h1 = line; h2 = line; h3 = line;

while t < 10000000

    delete(h0);
    delete(h1);
    delete(h2);
    delete(h3);
   
    hold on;
    h0 = line([x-r x+r],[y+r y+r],'Color','r');
    h1 = line([x-r x+r],[y-r y-r],'Color','r');
    h2 = line([x+r x+r],[y+r y-r],'Color','r');
    h3 = line([x-r x-r],[y+r y-r],'Color','r');

    was_a_key = waitforbuttonpress;
    
    if was_a_key && strcmp(get(fig, 'CurrentKey'), 'leftarrow')
        x = x - 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'rightarrow')
        x = x + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'a')
        r = r - 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'uparrow')
        y = y - 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'downarrow')
        y = y + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 's')
        r = r + 1;
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'return')
        FOV_coords = [x y r*2];
        break  
    elseif was_a_key && strcmp(get(fig, 'CurrentKey'), 'escape')
        close(fig);break
        
    end
end

%%
realWindowRadius = 2; %mm
scaleFactor = 2/coords(3);
real_window_center = [0 1];
FOV_center = [FOV_coords(1)-coords(1) -(FOV_coords(2)-coords(2))];
real_FOV_mm = [FOV_center*scaleFactor + real_window_center]


%%
figure;
[~, idx] = unique(coordsList(:,1));
uniqueFOVs = coordsList(idx,:);

winx = 0;
winy = 1;
winr = 2;

th = 0:pi/50:2*pi;
xunit = winr * cos(th) + winx;
yunit = winr * sin(th) + winy;
h = plot(xunit, yunit,'k','LineWidth',2);

r = 0.840/2;

for f = 1:length(uniqueFOVs)
    x = uniqueFOVs(f,1);
    y = uniqueFOVs(f,2);
    hold on;
    h1 = line([x-r x+r],[y+r y+r],'Color',[.8 1 0],'LineWidth',2);
    h2 = line([x-r x+r],[y-r y-r],'Color',[.8 1 0],'LineWidth',2);
    h3 = line([x+r x+r],[y+r y-r],'Color',[.8 1 0],'LineWidth',2);
    h4 = line([x-r x-r],[y+r y-r],'Color',[.8 1 0],'LineWidth',2);
end

axis square
box off
axis off
set(gca,'tickdir','out')