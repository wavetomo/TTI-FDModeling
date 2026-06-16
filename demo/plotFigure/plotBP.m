clc, clear, close all;
addpath('./lib');
%%
nz = 1126;
nx = 1251;
nt = 8001;
dx = 10;
dz = 10;
dt = 0.001;
x = (0:nx - 1) * dx / 1000;
z = (0:nz - 1) * dz / 1000;
t = (0:nt - 1) * dt;

%%
fontName = 'Arial';
fontSize = 9;
type = 3;

figNum = 9;

%%
vp = read2bin('../test_bp01/vp.bin', nz, nx);
epsilon = read2bin('../test_bp01/epsilon.bin', nz, nx);
delta = read2bin('../test_bp01/delta.bin', nz, nx);
theta = -read2bin('../test_bp01/theta.bin', nz, nx);

p1 = read2bin('../test_bp01/snap_p_01126_01251_03000ms.bin', nz, nx);
p2 = read2bin('../test_bp02/snap_p_01126_01251_03000ms.bin', nz, nx);
p3 = read2bin('../test_bp03/snap_p_01126_01251_03000ms.bin', nz, nx);

q1 = read2bin('../test_bp01/snap_p_01126_01251_04000ms.bin', nz, nx);
q2 = read2bin('../test_bp02/snap_p_01126_01251_04000ms.bin', nz, nx);
q3 = read2bin('../test_bp03/snap_p_01126_01251_04000ms.bin', nz, nx);

clip = max(abs([p2; p3; q2; q3]), [], 'all') * 0.5;

%%
width = 12;
height = 10;

figNum = figNum + 1;
figName = sprintf('Fig%02d', figNum);

fh = figure(figNum);
th = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
fh.Units = 'centimeters';
fh.Position(3) = width;
fh.Position(4) = height;

set(groot, 'defaultTextFontName', fontName);
set(groot, 'defaultAxesFontName', fontName);
set(groot, 'defaultTextFontSize', fontSize);
set(groot, 'defaultAxesFontSize', fontSize);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultAxesLineWidth', 1.0);


nexttile, imagesc(x, z, vp), axis image;
colormap(gca, jet), colorbar;
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

nexttile, imagesc(x, z, theta), axis image;
colormap(gca, jet), colorbar;
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');


nexttile, imagesc(x, z, epsilon), axis image;
colormap(gca, jet), colorbar;
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');


nexttile, imagesc(x, z, delta), axis image;
colormap(gca, jet), colorbar;
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

print(gcf, [figName, '.pdf'], '-dpdf', '-r600');
%%
width = 16;
height = 9;

figNum = figNum + 1;
figName = sprintf('Fig%02d', figNum);

figure(figNum);
th = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
fh.Units = 'centimeters';
fh.Position(3) = width;
fh.Position(4) = height;

set(groot, 'defaultTextFontName', fontName);
set(groot, 'defaultAxesFontName', fontName);
set(groot, 'defaultTextFontSize', fontSize);
set(groot, 'defaultAxesFontSize', fontSize);
set(groot, 'defaultLineLineWidth', 1.0);
set(groot, 'defaultAxesLineWidth', 1.0);

nexttile, imagesc(x, z, p3), axis image;
colormap(gca, gray), clim([-1, 1]*clip);
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

nexttile, imagesc(x, z, p2), axis image;
colormap(gca, gray), clim([-1, 1]*clip);
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

nexttile, imagesc(x, z, p1), axis image;
colormap(gca, gray), clim([-1, 1]*clip);
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

nexttile, imagesc(x, z, q3), axis image;
colormap(gca, gray), clim([-1, 1]*clip);
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

nexttile, imagesc(x, z, q2), axis image;
colormap(gca, gray), clim([-1, 1]*clip);
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

nexttile, imagesc(x, z, q1), axis image;
colormap(gca, gray), clim([-1, 1]*clip);
xlabel('Distance / km'), ylabel('Depth / km');
set(gca, 'XAxisLocation', 'top');
set(gca, 'Units', 'centimeters');
set(gca, 'TickDir', 'out');
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on');

cb = colorbar;
cb.Layout.Tile = 'east';

print(gcf, [figName, '.pdf'], '-dpdf', '-r600');

%%
