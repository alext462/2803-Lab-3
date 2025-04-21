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
function [x, t]=CSYSTEMlsim(Kptheta, Kdtheta, amplitude, overshoot_line)
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
    t=0:0.01:10;
    u=ones(1,length(t));
    [x,t] = lsim(sysTF, u, t);
x=x*amplitude; 

% 5 percent settling time
for j = 1:length(x)
    if x(j) >= 0.5*.05+0.5 || x(j) <= 0.5-0.05*0.5
        settling_time = t(j)
        x(j)
    end
end

figure(i)
grid on
plot(t, x)
titlename=['Lsim Set ' num2str(i)];
title(titlename)
xlabel('Time (s)')
ylabel('Theta (rad)')
PlotName = ['Lsim Set ' num2str(i)];
print(PlotName,'-dpng','-r300')
yline(overshoot_line,'r--', 'LineWidth', 1.5);
hold off;

max(x)
end 

end 
