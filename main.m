clc;
clear;
close all;

% Load the data
data1 = readmatrix('17and1');
data3 = readmatrix('7and1');

% Extract columns (assuming 2 columns: time [ms] and value)
time1 = data1(:,1); value1 = data1(:,2);
time3 = data3(:,1); value3 = data3(:,2);

% Normalize time to start at zero, convert to seconds, and shift by +2 sec
time1 = (time1 - min(time1)) / 1000 + 2;
time3 = (time3 - min(time3)) / 1000 + 2 + 1.1; % additional +1 second shift

% Keep only values where time <= 6 seconds
% mask1 = time1 <= 6;
% mask3 = time3 <= 6;

% time1 = time1(mask1);
% value1 = value1(mask1);
% 
% time3 = time3(mask3);
% value3 = value3(mask3);

%% Plot 1: K1 = 17, K3 = 1 for 2.4e
figure(1);
plot(time1, value1, 'b-', 'LineWidth', 1.5); hold on;
plot([0 5], [0.5 0.5], '--k'); % Reference line +0.5
plot([5 10], [-0.5 -0.5], '--k'); % Reference line -0.5
plot([5 5], [0.5 -0.5], '--k'); % Transition line
xlabel('Time (s)');
ylabel('Measured Value');
title('K_1 = 17, K_3 = 1');
% xlim([2, 6]);
grid on;
print('17 and 1 Hardware', '-dpng', '-r300')

%% Plot 2: K1 = 7, K3 = 1 for 2.4c
figure(2);
plot(time3, value3, 'g-', 'LineWidth', 1.5); hold on;
plot([0 5], [0.5 0.5], '--r'); % Reference line +0.5
plot([5 10], [-0.5 -0.5], '--b'); % Reference line -0.5
plot([5 5], [0.5 -0.5], '--k'); % Transition line
xlabel('Time (s)');
ylabel('Measured Value');
title('K_1 = 7, K_3 = 1');
% xlim([3, 6]); % Adjusted x-limits since data shifted by +3
grid on;
print('7 and 1 Hardware', '-dpng', '-r300')