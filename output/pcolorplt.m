clc; clear all

% d1 = importdata('inc_mag.txt');
% d2 = importdata('inc_phase.txt');

% d1 = importdata('oam1.txt');
% d2 = importdata('oam2.txt');

% d1 = importdata('pencil1.txt');
% d2 = importdata('pencil2.txt');

d1 = importdata('t30p0m1_mag.txt');
d2 = importdata('t30p0m1_phase.txt');

figure;
pcolor(d1);
shading flat;
colormap jet;
axis image;

figure;
pcolor(d2);
shading flat;
colormap jet;
axis image;