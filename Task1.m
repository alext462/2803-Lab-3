clc;
clear all;
close all;

% NACA = '2421';
% m = str2num(NACA(1));
% p = str2num(NACA(2));
% t = str2num(NACA(3:4));
m = 2;
p = 4;
t = 21;

n = 50;
c = 100;

%incr = 2;
%n = 180/incr;

theta = linspace(0, 180, n);
x = (c/2)*(1+cosd(theta));

%% ASEN 3802 - Lab 3 Part 1 - Main
% Provide a breif summary of the problem statement and code so that
% you or someone else can later understand what you attempted to do
% it doesn't have to be that long.
%
% Author: {Aleksei Suchkov}
% Collaborators: 
% Date: 

clc;
clear;
close all;

% Define parameters for the NACA 0021 airfoil
c = 100;  % Chord length
m = [0,2/c] ; % Maximum camber
p = [0,4/10];  % Location of maximum camber
t = [0.21, 0.21]; % Maximum thickness

N = 20;  % Number of points

% Generate the airfoil coordinates using the NACA_Airfoils function


for i = 1:2

[x_b, y_b] = NACA_Airfoils(m(i), p(i), t(i), c, N);
figure;
plot(x_b, y_b);
ylim([-40, 40]);
end

function [x_b,y_b] = NACA_Airfoils(m,p,t,c,N)


%incr = 2;
%n = 180/incr;

theta = linspace(0, 180, N);
x = (c/2)*(1+cosd(theta));

%x = linspace(0, c, N);
%disp(x)
y_t = ((t * c) / 0.2) * (0.2969 * sqrt(x / c) - 0.1260 * (x / c) - 0.3516 * (x / c).^2 + 0.2843 * (x / c).^3 - 0.1036 * (x / c).^4);

if all(x <= 0) || all(x < p*c)
    y_c = ((m * x) / p^2) * (2*p - (x / c));
    dyc_dx = (2*m) / p - (2*x*m) / (p^2 * c);
elseif all(x <= p*c) || all(x <= c)
    y_c = m * ((c - p) / (1 - p)^2) * (1 + (x / c) - 2*p);
    dyc_dx = -m / (1-p)^2 + m * (c - 2*x) / ((1-p)^2 * c) + (2*m*p / (1-p)^2);
end

zeta = atan(dyc_dx);

x_U = x - y_t.*sin(zeta);
x_L = x + y_t.*sin(zeta);

y_U = y_c + y_t.*cos(zeta);
y_L = y_c - y_t.*cos(zeta);

x_b = [fliplr(x_L), x_U];
y_b = [fliplr(y_L), y_U];

end
