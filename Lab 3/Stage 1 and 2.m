%% Week 6: LAB 3 4
% Name: Avvienash A/L Jaganathan
% ID: 32281013
% Date: 8/4/2022
clear all; close all; clc;

% constants

A =1;
lamda = 1;
m = 1;
g = 9.81;
Ri = 3;
dt = 2e-4;
k = 2*pi/lamda;

n =20000;
% initialise
t = 0:dt:dt*(n-1);
q = zeros(3,n);
q(2,1) = 7;


for a = 2:n
    % use eq motion to find ac
    term1 = m*(1+(A^2)*(k^2)*(sin(k*q(1,a-1)))^2);
    term2 = m*q(2,a-1)^2 * A^2 * k^3 * sin(k*q(1,a-1))* cos(k*q(1,a-1));
    term3 = m*g*A*k*sin(k*q(1,a-1));
    term4 = Ri*q(2,a-1);
    q(3,a)= -(term2+term3+term4)/term1;
    
    % use a to calc rest
    q(2,a) = q(2,a-1) + q(3,a)*dt;
    q(1,a) = q(1,a-1)+q(2,a-1)*dt + q(3,a)*(dt^2)*0.5;
end

plot(t,q(1,:),t,q(2,:),t,q(3,:));
legend("q","qd","qdd");

figure(2)
plot(q(1,:),q(2,:));
xlabel("q")
ylabel("qd")
title("phase")

x = q(1,:);
vx = q(2,:);
vy = A.*k.*vx.*sin(k.*x);

figure(3)
plot(x,sqrt(vx.^2 + vy.^2));
xlabel("q")
ylabel("abs V")
title("Absolute Velocity")