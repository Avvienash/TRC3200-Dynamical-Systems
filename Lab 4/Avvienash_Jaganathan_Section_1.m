%% Week 6: LAB 4
% Name: Avvienash A/L Jaganathan
% ID: 32281013
% Date: 29/4/2022
clear all; close all; clc;

%% Stage 1
% constants


m = 1;
l = 3;
g = 9.81;
k = 2;
R = 0.5;
dt = 2e-4;

n =20000;
% initialise
t = 0:dt:dt*(n-1);
q = zeros(3,n);
q(2,1) = 8;


for a = 2:n

    % use eq motion to find ac
    term1 = m*l*l;
    term2 = m*g*l*sin(q(1,a-1));
    term3 = 4*k*(R^2)*sin(q(1,a-1))*cos(q(1,a-1));
    q(3,a)= -(term2+term3)/term1;
    
    % use a to calc rest
    q(2,a) = q(2,a-1) + q(3,a)*dt;
    q(1,a) = (q(1,a-1)+q(2,a-1)*dt + q(3,a)*(dt^2)*0.5);
end

figure(1)
plot(q(1,:),q(2,:));
xlabel("theta")
ylabel("d.theta")
title("Stage 1: Phase-Plane");

%% Stage 2

rf = 2;

d = zeros(3,n);
d(2,1) = 8;
V2 = zeros(1,n);

for a = 2:n

    % use eq motion to find ac
    term1 = m*l*l;
    term2 = m*g*l*sin(d(1,a-1));
    term3 = 4*k*(R^2)*sin(d(1,a-1))*cos(d(1,a-1));
    term4 = rf*d(2,a-1);
    d(3,a)= -(term2+term3+term4)/term1;
    
    % use a to calc rest
    d(2,a) = d(2,a-1) + d(3,a)*dt;
    d(1,a) = (d(1,a-1)+d(2,a-1)*dt + d(3,a)*(dt^2)*0.5);
    V2(1,a) = m*g*l*(1-cos(d(1,a))) + 2*k*R*R*(sin(d(1,a)) )^2;
end

figure(2)
plot(d(1,:),d(2,:),q(1,:),q(2,:));
legend("Rf = 2","Rf = 0");
xlabel("Theta")
ylabel("d.theta")
title("Stage 2: phase-plane");

%% Stage 3

d3 = zeros(3,n);
V3 = zeros(1,n);
d3(2,1) = 8;


for a = 2:n

    % use eq motion to find ac
    term1 = m*l*l;
    term2 = m*g*l*sin(d3(1,a-1));
    term3 = -2*k*(R^2)*sin(2*(d3(1,a-1) - (pi/2)));
    term4 = 4*k*(R^2)*sin(d3(1,a-1) - (pi/2));
    term5 = rf*d3(2,a-1);
    d3(3,a)= -(term2+term3+term4+term5)/term1;
    
    % use a to calc rest
    d3(2,a) = d3(2,a-1) + d3(3,a)*dt;
    d3(1,a) = (d3(1,a-1)+d3(2,a-1)*dt + d3(3,a)*(dt^2)*0.5);
    V3(1,a) = m*g*l*(1-cos(d3(1,a))) + 2*k*R*R*(1 + sin(d3(1,a)) )^2;

end

figure(3)
plot(d3(1,:),V3(1,:),d(1,:),V2(1,:));
legend("Stage 3","Stage 2");
xlabel("Theta")
ylabel("Potential Energy,V")
title("Stage 3: Potential Eneregy agaisnt Theta");


R = R*2;

q32 = zeros(3,n);
q32(2,1) = 8;
V32 = zeros(1,n);

for a = 2:n

    % use eq motion to find ac
    term1 = m*l*l;
    term2 = m*g*l*sin(q32(1,a-1));
    term3 = 4*k*(R^2)*sin(q32(1,a-1))*cos(q32(1,a-1));
    term4 = rf*q32(2,a-1);
    q32(3,a)= -(term2+term3+term4)/term1;
    
    % use a to calc rest
    q32(2,a) = q32(2,a-1) + q32(3,a)*dt;
    q32(1,a) = (q32(1,a-1)+q32(2,a-1)*dt + q32(3,a)*(dt^2)*0.5);
    V32(1,a) = m*g*l*(1-cos(q32(1,a))) + 2*k*R*R*(sin(q32(1,a)) )^2;
end

q33 = zeros(3,n);
V33 = zeros(1,n);
q33(2,1) = 8;


for a = 2:n

    % use eq motion to find ac
    term1 = m*l*l;
    term2 = m*g*l*sin(q33(1,a-1));
    term3 = -2*k*(R^2)*sin(2*(q33(1,a-1) - (pi/2)));
    term4 = 4*k*(R^2)*sin(q33(1,a-1) - (pi/2));
    term5 = rf*q33(2,a-1);
    q33(3,a)= -(term2+term3+term4+term5)/term1;
    
    % use a to calc rest
    q33(2,a) = q33(2,a-1) + q33(3,a)*dt;
    q33(1,a) = (q33(1,a-1)+q33(2,a-1)*dt + q33(3,a)*(dt^2)*0.5);
    V33(1,a) = m*g*l*(1-cos(q33(1,a))) + 2*k*R*R*(1 + sin(q33(1,a)) )^2;

end

figure(4)
plot(q33(1,:),V33(1,:),q32(1,:),V32(1,:));
legend("Stage 3","stage 2");
xlabel("theta")
ylabel("V")
title("Potential Energy (Double R) ");