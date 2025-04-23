clc
clear
close all

% Preallocate vectors

    amplitude=0.5;

    overshoot_line = amplitude * 0.2+amplitude;
    

    Kptheta = [7];
    
    Kdtheta = [1];
   
    
%[x, t]=CSYSTEMstep(Kptheta, Kdtheta, amplitude);

[x2, t2]=CSYSTEMlsim(Kptheta, Kdtheta, amplitude, overshoot_line);

%% Functions 

clc
clear
close all

% Preallocate vectors

    amplitude=0.5; 
    

    Kptheta = [10 20 5 10 10 10 ];
    
    Kdtheta = [0 0 0 1 -1 -0.5];
   
    
%[x, t]=CSYSTEMstep(Kptheta, Kdtheta, amplitude);

[x2, t2]=CSYSTEMlsim(Kptheta, Kdtheta, amplitude)

%% Functions 

% 1
function [x, t]=CSYSTEMstep(Kptheta, Kdtheta, amplitude)
Kg = 33.3;
Km = 0.0401;
J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
Rm = 19.2;

for i=1:length(Kdtheta)
n1 = Kptheta(i) * Kg * Km / (J * Rm);
    d2 = 1;
    d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(i)* Kg* Km / (J * Rm);
    d0 = Kptheta(i) * Kg * Km / (J * Rm);
% Define transfer function
    num = n1;
    den = [d2 d1 d0];
    sysTF = tf(num, den);
    % Step Response
    [x1,t1] = step(sysTF);
    [x2,t2]=step(sysTF);
x1=x1*amplitude; 
x2=-x2*amplitude;
figure(i)
grid on
x=[x1:x2];
t=[t1:t2];
plot(t, x)
titlename=['Step Set ' num2str(i)];
title(titlename)
xlabel('Time (s)')
ylabel('Theta (rad)')
PlotName = ['Step Set ' num2str(i)];
print(PlotName,'-dpng','-r300')
end 

end 

% 2
% function [x, t]=CSYSTEMlsim(Kptheta, Kdtheta, amplitude)
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
%     t1=0:0.01:5;
%     t2=5:0.01:10;
%     u=ones(1,length(t1));
%     [x1,t_1] = lsim(sysTF, u, t1);
%     [x2, t_2]=lsim(sysTF, u, t2)
% x1=x1*amplitude; 
% x2=-x2*amplitude;
% 
% figure(i)
% grid on
% plot(t, x)
% titlename=['Lsim Set ' num2str(i)];
% title(titlename)
% xlabel('Time (s)')
% ylabel('Theta (rad)')
% PlotName = ['Lsim Set ' num2str(i)];
% print(PlotName,'-dpng','-r300')
% end 
% 
% end 

% function [x, t] = CSYSTEMlsim(Kptheta, Kdtheta, amplitude)
% Kg = 33.3;
% Km = 0.0401;
% J = 0.0005 + 0.2 * (0.2794)^2 + 0.0015;
% Rm = 19.2;
% 
% t1 = 0:0.01:5;
% t2 = 5.01:0.01:10;
% t = [t1 t2];
% u1 = amplitude * ones(size(t1));
% u2 = amplitude * ones(size(t2));
% 
% for i = 1:length(Kdtheta)
%     n1 = Kptheta(i) * Kg * Km / (J * Rm);
%     d2 = 1;
%     d1 = Kg^2 * Km^2 / (J * Rm) + Kdtheta(i) * Kg * Km / (J * Rm);
%     d0 = Kptheta(i) * Kg * Km / (J * Rm);
% 
%     % Define transfer function
%     num = n1;
%     den = [d2 d1 d0];
%     sysTF = tf(num, den);
% 
%     % Simulate forward response
%     x1 = lsim(sysTF, u1, t1);
%     % Simulate reverse response
%     x2 = -lsim(sysTF, u2, t2);
% 
%     % Combine responses
%     x = [x1; x2];
% 
%     % Plot
%     figure(i)
%     plot(t, x)
%     grid on
%     title(['Lsim Set ' num2str(i)])
%     xlabel('Time (s)')
%     ylabel('Theta (rad)')
% 
%     % Save plot
%     PlotName = ['Lsim Set ' num2str(i)];
%     print(PlotName, '-dpng', '-r300')
% end
% end
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
    title(['Lsim Set ' num2str(i)])
    xlabel('Time (s)')
    ylabel('Theta (rad)')
% Horizontal lines
plot([0 5], [0.5 0.5], '--r', 'DisplayName', 'Ref +0.5')
plot([5 10], [-0.5 -0.5], '--b', 'DisplayName', 'Ref -0.5')

% Vertical connector
plot([5 5], [0.5 -0.5], '--k', 'DisplayName', 'Transition')

    % Save plot
    PlotName = ['Lsim Set ' num2str(i)];
    print(PlotName, '-dpng', '-r300')
end
end

% 1
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
%     [x,t] = step(sysTF);
% x=x*amplitude; 
% figure(i)
% grid on
% plot(t, x)
% titlename=['Step Set ' num2str(i)];
% title(titlename)
% xlabel('Time (s)')
% ylabel('Theta (rad)')
% PlotName = ['Step Set ' num2str(i)];
% print(PlotName,'-dpng','-r300')
% end 
% 
% end 

% 2
% function [x, t]=CSYSTEMlsim(Kptheta, Kdtheta, amplitude, overshoot_line)
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
%     t=0:0.01:10;
%     u=ones(1,length(t));
%     [x,t] = lsim(sysTF, u, t);
% x=x*amplitude; 
% 
% % 5 percent settling time
% for j = 1:length(x)
%     if x(j) >= 0.5*.05+0.5 || x(j) <= 0.5-0.05*0.5
%         settling_time = t(j)
%         x(j)
%     end
% end
% 
% figure(i)
% grid on
% plot(t, x)
% titlename=['Lsim Set ' num2str(i)];
% title(titlename)
% xlabel('Time (s)')
% ylabel('Theta (rad)')
% PlotName = ['Lsim Set ' num2str(i)];
% print(PlotName,'-dpng','-r300')
% yline(overshoot_line,'r--', 'LineWidth', 1.5);
% hold off;
% 
% max(x)
% end 
% 
% end 
