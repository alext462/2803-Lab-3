clc
clear
close all

% Preallocate vectors

    amplitude=0.5; 
    
    Kptheta = [7, 17];
    
    Kdtheta = [1, 1];

%% Functions 

data1 = readmatrix('17and1');
data3 = readmatrix('7and1');

% Extract columns (assuming 2 columns: time [ms] and value)
time1 = data1(:,1); value1 = data1(:,2);
time3 = data3(:,1); value3 = data3(:,2);

% Normalize time to start at zero, convert to seconds, and shift by +2 sec
time1 = (time1 - min(time1)) / 1000+2;
time3 = (time3 - min(time3)) / 1000; % additional +1 second shift

t = 0:0.01:10; % Full time vector
u = amplitude * ones(size(t)); 
%% plot 1 for 2.4b


Kg = 33.3;  
Km = 0.0401;
J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
Rm = 19.2;


u(t >= 5) = -amplitude; % Step down at t = 5 to create a smooth transition


    n1 = Kptheta(1) * Kg * Km / (J * Rm);
    d2 = 1;
    d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(1) * Kg * Km / (J * Rm);
    d0 = Kptheta(1) * Kg * Km / (J * Rm);

    % Define transfer function
    num = n1;
    den = [d2 d1 d0];
    sysTF = tf(num, den);

    % Simulate system response
    x = lsim(sysTF, u, t);

    % Plot
    figure(1)
    hold on
    plot(t, x)
    % plot(time1, value1, 'b-', 'LineWidth', 1.5); 
    
    grid on
    title(['K_p = 7 , K_d = 1'])
    xlabel('Time (s)')
    ylabel('Theta (rad)')
% Horizontal lines
plot([0 5], [0.5 0.5], '--k', 'DisplayName', 'Ref +0.5')
plot([5 10], [-0.5 -0.5], '--k', 'DisplayName', 'Ref -0.5')

% Vertical connector
plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')

    % Save plot
   
    print('K_p = 7 , K_d = 1 Model', '-dpng', '-r300')
hold off


%% plot 2



u(t >= 5) = -amplitude; % Step down at t = 5 to create a smooth transition

%for i = 1:length(Kdtheta)
    n1 = Kptheta(2) * Kg * Km / (J * Rm);
    d2 = 1;
    d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(2) * Kg * Km / (J * Rm);
    d0 = Kptheta(2) * Kg * Km / (J * Rm);

    % Define transfer function
    num = n1;
    den = [d2 d1 d0];
    sysTF = tf(num, den);

    % Simulate system response
    x = lsim(sysTF, u, t);

    % Plot
    figure(2)
    hold on
    plot(t, x)
    % plot(time3, value3, 'g-', 'LineWidth', 1.5); hold on;
    grid on
    title(['K_p = 17 , K_d = 1'])
    xlabel('Time (s)')
    ylabel('Theta (rad)')
% Horizontal lines
plot([0 5], [0.5 0.5], '--k', 'DisplayName', 'Ref +0.5')
plot([5 10], [-0.5 -0.5], '--k', 'DisplayName', 'Ref -0.5')

% Vertical connector
plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')

    % Save plot
    
    print('K_p = 17 , K_d = 1 Model', '-dpng', '-r300')
hold off

%% Plot 3 for 3.1
% Kg = 33.3;  
% Km = 0.0401;
% J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
% Rm = 19.2;
% 
% 
% u(t >= 5) = amplitude; % Step down at t = 5 to create a smooth transition
% 
%     n1 = Kptheta(2) * Kg * Km / (J * Rm);
%     d2 = 1;
%     d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(2) * Kg * Km / (J * Rm);
%     d0 = Kptheta(2) * Kg * Km / (J * Rm);
% 
%     % Define transfer function
%     num = n1;
%     den = [d2 d1 d0];
%     sysTF = tf(num, den);
% 
%     % Simulate system response
%     x = lsim(sysTF, u, t);
% 
%     % algorithm to find 5% settling time
% 
%     upper=0.5+0.05*0.5;
%     lower=0.5-0.05*0.5;
% 
%     for i = 1: length(t)
%         if x(i) > upper || x(i)<lower && t(i)<4
%             settlingtime=t(i); 
%             settlingx=x(i);
%         end 
%   if x(i)==max(x)
%       max_x=x(i);
%       max_t=t(i);
%   end 
% 
%     end 
% 
% 
% 
%     % Plot
%     figure(3)
%     hold on
%     xline(settlingtime, '--r', 'LineWidth', 1.5)
%     yline(max_x, '--p', 'LineWidth', 1.5)
%     plot(t, x)
%     plot(time1, value1, 'b-', 'LineWidth', 1.5); 
% 
%     grid on
%     title(['K_p = 17 , K_d = 1 Comparison'])
%     xlabel('Time (s)')
%     ylabel('Theta (rad)')
% % Horizontal lines
% plot([0 5], [0.5 0.5], '--k', 'DisplayName', 'Ref +0.5')
% plot([5 10], [-0.5 -0.5], '--k', 'DisplayName', 'Ref -0.5')
% 
% % Vertical connector
% plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')
% 
%     % Save plot
% 
%     print('K_p = 17 , K_d = 1 Comparison', '-dpng', '-r300')
% hold off
% 


Kg = 33.3;  
Km = 0.0401;
J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
Rm = 19.2;

u(t >= 5) = amplitude; % Step down at t = 5 to create a smooth transition

n1 = Kptheta(2) * Kg * Km / (J * Rm);
d2 = 1;
d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(2) * Kg * Km / (J * Rm);
d0 = Kptheta(2) * Kg * Km / (J * Rm);

% Define transfer function
num = n1;
den = [d2 d1 d0];
sysTF = tf(num, den);

% Simulate system response
x = lsim(sysTF, u, t);

% algorithm to find 5% settling time of model
upper = 0.5 + 0.05 * 0.5;
lower = 0.5 - 0.05 * 0.5;

for i = 1:length(t)
    if (x(i) > upper || x(i) < lower) 
        settlingtimeModel = t(i); 
        settlingx = x(i);
    end 
    if x(i) == max(x)
        max_x = x(i);
        max_t = t(i);
    end 
end

% algorithm to find 5% settling time of hardware
counter=0; 
for j = 1:length(time1)
    if (value1(j) > upper || value1(j) < lower) 
        settlingtimeHardware = time1(j); 
        settlingValue = value1(j);
    end 
    if value1(j) == max(value1) && counter==0
        max_Value = value1(j);
        max_tvalue = time1(j);
        counter=counter+1;
    end 
end

% Plot
figure(3)
hold on
plot(time1, value1, 'b-', 'LineWidth', 1.5); 
plot(t, x, 'LineWidth', 1.5)
 xline(settlingtimeModel, '--r', 'LineWidth', 1.5)
 xline(settlingtimeHardware,'--b', 'LineWidth', 1.5)
plot(max_tvalue, max_Value, '*')
plot(max_t, max_x, '*')

% Label overshoot and 5% settling time
% text(max_t, max_x, ['  Overshoot: ', num2str(max_x, '%.2f')], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold')
% text(settlingtime, settlingx, ['  5% Settling at t = ', num2str(settlingtime, '%.2f') ' s'], 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 10, 'FontWeight', 'bold')
% Label overshoot and 5% settling time with coordinate offsets
text(max_t + 0.3, max_x-0.01, ['Overshoot: ', num2str(max_x, '%.2f')], ...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left', ...
    'FontSize', 7)

text(settlingtimeModel + 0.2, settlingx - 0.05, ['5% Settling at t = ', num2str(settlingtimeModel, '%.2f'), ' s'], ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    'FontSize', 7)

% Horizontal reference lines
plot([0 10], [0.5 0.5], '--k', 'DisplayName', 'Ref +0.5')
% plot([5 10], [-0.5 -0.5], '--k', 'DisplayName', 'Ref -0.5')
% plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')

grid on
title('K_p = 17 , K_d = 1 Comparison')
xlabel('Time (s)')
ylabel('Theta (rad)')

print('K_p = 17 , K_d = 1 Comparison', '-dpng', '-r300')

hold off


%% plot 4 for 2.4d

u(t >= 5) = -amplitude; % Step down at t = 5 to create a smooth transition


    n1 = Kptheta(1) * Kg * Km / (J * Rm);
    d2 = 1;
    d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(1) * Kg * Km / (J * Rm);
    d0 = Kptheta(1) * Kg * Km / (J * Rm);

    % Define transfer function
    num = n1;
    den = [d2 d1 d0];
    sysTF = tf(num, den);

    % Simulate system response
    x = lsim(sysTF, u, t);

    % Plot
    figure(4)
    hold on
    plot(t, x)
     plot(time3, value3, 'g-', 'LineWidth', 1.5); hold on;
    grid on
    title(['K_p = 7 , K_d = 1 Comparison'])
    xlabel('Time (s)')
    ylabel('Theta (rad)')
% Horizontal lines
plot([0 5], [0.5 0.5], '--k', 'DisplayName', 'Ref +0.5')
plot([5 10], [-0.5 -0.5], '--k', 'DisplayName', 'Ref -0.5')

% Vertical connector
plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')

    % Save plot
    
    print('K_p = 7 , K_d = 1 Comparison', '-dpng', '-r300')
hold off
