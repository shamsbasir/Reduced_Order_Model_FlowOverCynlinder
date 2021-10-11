function f1 = plotCylinder(VORT,C)

set(0,'defaulttextinterpreter','latex')
% f1 = figure
vortmin = -5;  % only plot what is in -5 to 5 range
vortmax = 5;
VORT(VORT>vortmax) = vortmax;  % cutoff at vortmax
VORT(VORT<vortmin) = vortmin;  % cutoff at vortmin

imagesc(VORT); % plot vorticity field
load CCcool.mat 
colormap(CC);  % use custom colormap

% clean up axes
% set(gca,'XTick',[1 50 100 150 200 250 300 350 400 449],'XTickLabel',{'-1','0','1','2','3','4','5','6','7','8'})
% set(gca,'YTick',[1 50 100 150 199],'YTickLabel',{'2','1','0','-1','-2'});
% set(gcf,'Position',[500 500 600 260])

set(gca,'XTick',[1 50 100 150 200 250 300 350 400 450 500],'XTickLabel',{'-2','-1','0','1','2','3','4','5','6','7','8'})
set(gca,'YTick',[1 50 100 150 200 250],'YTickLabel',{'-2.5','-1.5','-.5','.5','1.5','2.5'});
set(gcf,'Position',[500 500 1000 520])

axis equal
hold on

% add contour lines (positive = solid, negative = dotted)
contour(VORT,C,'.k','LineWidth',1.2)
% contour(VORT,[.125 .25 .5:.5:5.5],'-k','LineWidth',1.2)

theta = (1:100)/100'*2*pi;
% x = 49+25*sin(theta);
% y = 99+25*cos(theta);
x = 101+25*sin(theta);
y = 126+25*cos(theta);
fill(x,y,[.3 .3 .3])  % place cylinder
plot(x,y,'k','LineWidth',1.2) % cylinder boundary

set(gcf,'PaperPositionMode','auto') % 