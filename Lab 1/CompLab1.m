%% Week 2: LAB 1
% Name: Avvienash A/L Jaganathan
% ID: 32281013
% Date: 6/3/2022
clear all; close all; clc;
%% Stage 1
n = 2e7;
a = zeros(1,n);
for loop = 1:n
    a(loop) = rand;
end
average = sum(a)/n;
disp(['The average is ' num2str(average)]) ;
%% Stage 2
x=zeros(1,n);
y=zeros(1,n);
z=zeros(1,n);
inside=zeros(1,n);

for loop=1:n
    x(loop)=2*rand-1;
    y(loop)=2*rand-1;
    z(loop)=2*rand-1;
    
    if sqrt(x(loop)^2+y(loop)^2 +z(loop)^2)<=1
        inside(loop)=1;
    else
        inside(loop)=0;
    end
end
%% Stage 3
volest = sum(inside)*8/n;
disp(['The estimated volume is ' num2str(volest)]) 
%% Stage 4
V = 4/3 * pi;
disp(['The actual volume is ' num2str(V)])
Error = abs((V-volest)/V *100);
disp(['The error is ' num2str(Error) ' percent']) 
%% Stage 5
N = 1:1000;
VolEst = zeros(1,length(N));
for loop = 1:length(N)
    n =N(loop);
    x=zeros(1,n);
    y=zeros(1,n);
    z=zeros(1,n);
    inside=zeros(1,n);

    for loop=1:n
        x(loop)=2*rand-1;
        y(loop)=2*rand-1;
        z(loop)=2*rand-1;

        if sqrt(x(loop)^2+y(loop)^2 +z(loop)^2)<=1
            inside(loop)=1;
        else
            inside(loop)=0;
        end
    end
    VolEst(loop) = sum(inside)*8/n;
end
plot(N,VolEst)
