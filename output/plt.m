clc; clear all;

% d1 = importdata('t0p0.txt');
% d2 = importdata('t0p0_simu.csv');

% d1 = importdata('t0p0_90.txt');
% d2 = importdata('t0p0_90_simu.csv');

d1 = importdata('t30p0m1.txt');
d2 = importdata('t30p0m1_simu.csv');

d2 = d2.data;

figure;
plot(rad2deg(d1(:,1)), d1(:,2), d2(:,1), d2(:,2), ...
     'LineWidth', 2);
ylim([-20, 30]);
xlabel('Theta (Deg.)');
ylabel('Gainn (dB)');
legend('Theoretical', 'Simulation');

