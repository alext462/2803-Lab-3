clc
clear
close all

% Preallocate vectors

    amplitude=0.5; 
    

    Kptheta = [7, 17];
    
    Kdtheta = [1, 1];
   
    
%[x, t]=CSYSTEMstep(Kptheta, Kdtheta, amplitude);

%[x2, t2]=CSYSTEMlsim(Kptheta, Kdtheta, amplitude)

%% Functions 




% % 1
% function [x, t]=CSYSTEMstep(Kptheta, Kdtheta, amplitude)
% Kg = 33.3;
% Km = 0.0401;
% J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
% Rm = 19.2;
% 
% for i=1:length(Kdtheta)
% n1 = Kptheta(i) * Kg * Km / (J * Rm);
%     d2 = 1;
%     d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(i)* Kg* Km / (J * Rm);
%     d0 = Kptheta(i) * Kg * Km / (J * Rm);
% % Define transfer function
%     num = n1;
%     den = [d2 d1 d0];
%     sysTF = tf(num, den);
%     % Step Response
%     [x1,t1] = step(sysTF);
%     [x2,t2]=step(sysTF);
% x1=x1*amplitude; 
% x2=-x2*amplitude;
% figure(i)
% grid on
% x=[x1:x2];
% t=[t1:t2];
% plot(t, x)
% titlename=['Step Set ' num2str(i)];
% title(titlename)
% xlabel('Time (s)')
% ylabel('Theta (rad)')
% PlotName = ['Step Set ' num2str(i)];
% print(PlotName,'-dpng','-r300')
% end 

% end 

%function [x, t] = CSYSTEMlsim(Kptheta, Kdtheta, amplitude)
Kg = 33.3;
Km = 0.0401;
J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
Rm = 19.2;

t = 0:0.01:10; % Full time vector
u = amplitude * ones(size(t)); 
u(t >= 5) = -amplitude; % Step down at t = 5 to create a smooth transition

%for i = 1:length(Kdtheta)
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
    plot(t, x)
    hold on
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
    PlotName = [num2str(1)];
    print(PlotName, '-dpng', '-r300')
hold off
%function [x, t] = CSYSTEMlsim(Kptheta, Kdtheta, amplitude)
Kg = 33.3;
Km = 0.0401;
J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
Rm = 19.2;

t = 0:0.01:10; % Full time vector
u = amplitude * ones(size(t)); 
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
    plot(t, x)
    hold on
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
    PlotName = [num2str(2)];
    print(PlotName, '-dpng', '-r300')
hold off


function [x, t] = CSYSTEMlsim(Kptheta, Kdtheta, amplitude)
Kg = 33.3;
Km = 0.0401;
J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
Rm = 19.2;

t = 0:0.01:10; % Full time vector
u = amplitude * ones(size(t)); 
u(t >= 5) = -amplitude; % Step down at t = 5 to create a smooth transition

for i = 1:length(Kdtheta)
    n1 = Kptheta(i) * Kg * Km / (J * Rm);
    d2 = 1;
    d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(i) * Kg * Km / (J * Rm);
    d0 = Kptheta(i) * Kg * Km / (J * Rm);

    % Define transfer function
    num = n1;
    den = [d2 d1 d0];
    sysTF = tf(num, den);

    % Simulate system response
    x = lsim(sysTF, u, t);

    % Plot
    figure(i)
    plot(t, x)
    hold on
    grid on
    title(['Angular Position of Arm vs. Time '])
    xlabel('Time (s)')
    ylabel('Theta (rad)')
% Horizontal lines
plot([0 5], [0.5 0.5], '--k', 'DisplayName', 'Ref +0.5')
plot([5 10], [-0.5 -0.5], '--k', 'DisplayName', 'Ref -0.5')

% Vertical connector
plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')

    % Save plot
    PlotName = [num2str(i)];
    print(PlotName, '-dpng', '-r300')
end
end