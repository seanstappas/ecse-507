f = @(x,y) (1-x).^2 + 100*(y-x.^2).^2;
x = linspace(-1.5,1.1); y = linspace(-1.5,1.1);
[xx,yy] = meshgrid(x,y); ff = f(xx,yy);
levels = 5:5:600;
LW = 'linewidth'; FS = 'fontsize'; MS = 'markersize';
figure, contour(x,y,ff,levels,LW,1.2), colorbar
axis([-1.5 1.1 -0.1 1.1]), axis square, hold on