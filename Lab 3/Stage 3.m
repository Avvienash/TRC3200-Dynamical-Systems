%% Week 6:PRE LAB 3 4
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

n =40000;
% initialise
t = 0:dt:dt*(n-1);
q = zeros(3,n);
inV = 4:0.02:(12-0.02);
EndWell = zeros(1,400);
for b = 1:400
    q(2,1) = inV(b);


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
        
        if ((q(1,a) - q(1,a-1)) < 0)
            break
        end
    end
    EndWell(b) = q(1,a);
end
EndWell2 = floor(EndWell./lamda);

FinalAns = zeros(3,400);
FinalAns(1,:) = inV;
FinalAns(2,:) = EndWell;
FinalAns(3,:) = EndWell2;

plot((1:400),EndWell2);
xlabel("Initial Velocity");
ylabel("End Well");
title("End Wells")
