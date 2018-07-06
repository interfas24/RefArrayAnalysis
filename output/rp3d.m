clc; clear all;

name = 't0p0_gain3d';
% name = 't30p0m1_gain3d';

A = importdata([name, '.txt']);

N = 180;
th = -30;
for i = 1:N+1
    for j = 1:N+1
        if A(i, j) < th
            A(i, j) = th;
        end
    end
end

minA = min(min(A));
maxA = max(max(A));

A = A - minA;

t = -pi/2:pi./N:pi/2;
tp = t';
p = 0:2*pi./N:2*pi;
[T, P] = meshgrid(t, p);
X = A.*sin(T).*cos(P);
Y = A.*sin(T).*sin(P);
Z = A.*cos(T);

minX = min(min(X));
minY = min(min(Y));
maxX = max(max(X));
maxY = max(max(Y));

figure
x0=5;
y0=5;
width=4.8;
height=4;
set(gcf,'units','inches','position',[x0,y0,width,height]);

AxesH = axes;
drawnow;
InSet = get(AxesH, 'TightInset');
set(AxesH, 'Position', [InSet(1:2), 1-InSet(1)-InSet(3), 1-InSet(2)-InSet(4)])

surf(X,Y,Z);
view(150, 36);
shading interp
axis image;
colormap jet
% grid off;
axis([minX*1.2 maxX*1.2 minY*1.2 maxY*1.2 0 maxA-minA]);
xlabel('X'), ylabel('Y'), zlabel('Z');

set(gca, 'Color', 'none')


