%% Week 10: COMP LAB TEST
% Name: Avvienash A/L Jaganathan
% ID: 32281013
% Date: 13/5/2022
clear all; close all; clc;

%% Stage 1

%constants
m = 3;
l = 3;
k1 = 2;
k2 = 0;
g = 9.81;
dt = 2e-4;
tsteps = 40001;

% create Variables
x = zeros(3,tsteps);
theta = zeros(3,tsteps);
Time = zeros(1,tsteps);

% initial Values
theta(2,1) = 5;
x(2,1) = 5;

%loop in time
for a=2:tsteps
    %
    
    x(3,a) = (1/(m+m-m*((cos(theta(1,a-1)))^2)))* ( ( -m*g* sin(theta(1,a-1)) * cos(theta(1,a-1)) ) + ( m*cos(theta(1,a-1))*k2*theta(1,a-1) ) + (m*l*(theta(2,a-1))^2 *sin(theta(1,a-1))) + ( -k1*x(1,a-1) ));
    theta(3,a) = (-x(3,a)*cos(theta(1,a-1))/l) + ( sin(theta(1,a-1))*g/l) - (k2*theta(1,a-1)/l);

    x(2,a) = x(2,a-1) + x(3,a)*dt;
    x(1,a) = (x(1,a-1)+x(2,a-1)*dt + x(3,a)*(dt^2)*0.5);

    theta(2,a) = theta(2,a-1) + theta(3,a)*dt;
    theta(1,a) = (theta(1,a-1)+theta(2,a-1)*dt + theta(3,a)*(dt^2)*0.5);

    Time(a)=Time(a-1)+dt;
end

% plot vel
 figure(1)
 plot(Time(1,:) ,x(2,:),Time(1,:) ,theta(2,:));
 title("Velocity againt time Stage 1");
 legend("xdt","thetadt");
 xlabel("Time(s)");
ylabel("Velocity");

%% Stage 2

%constants
dt = 2e-4;
tsteps = 50001;

% recreate Variables
x = zeros(3,tsteps);
theta = zeros(3,tsteps);
Time = zeros(1,tsteps);

% initial Values
x(1,1) = 4;

%loop in time
for a=2:tsteps
    %
    
    x(3,a) = (1/(m+m-m*((cos(theta(1,a-1)))^2)))* ( ( -m*g* sin(theta(1,a-1)) * cos(theta(1,a-1)) ) + ( m*cos(theta(1,a-1))*k2*theta(1,a-1) ) + (m*l*(theta(2,a-1))^2 *sin(theta(1,a-1))) + ( -k1*x(1,a-1) ));
    theta(3,a) = (-x(3,a)*cos(theta(1,a-1))/l) + ( sin(theta(1,a-1))*g/l) - (k2*theta(1,a-1)/l);

    x(2,a) = x(2,a-1) + x(3,a)*dt;
    x(1,a) = (x(1,a-1)+x(2,a-1)*dt + x(3,a)*(dt^2)*0.5);

    theta(2,a) = theta(2,a-1) + theta(3,a)*dt;
    theta(1,a) = (theta(1,a-1)+theta(2,a-1)*dt + theta(3,a)*(dt^2)*0.5);

    Time(a)=Time(a-1)+dt;
end

% plot vel
 figure(2)
 plot(Time(1,:) ,x(2,:),Time(1,:) ,theta(2,:));
 title("Velocity againt time stage 2");
 legend("xdt","thetadt");
 xlabel("Time(s)");
 ylabel("Velocity");

 %% Stage 3

 % recreate Variables
x = zeros(3,tsteps);
theta = zeros(3,tsteps);
Time = zeros(1,tsteps);
tsteps = 100001;


% initial Values
x(1,1) = 4;
k2 = 2;

boolean = 1;
while (boolean == 1)
    k2 = k2 + 2;
    for a=2:tsteps
        %
        
        x(3,a) = (1/(m+m-m*((cos(theta(1,a-1)))^2)))* ( ( -m*g* sin(theta(1,a-1)) * cos(theta(1,a-1)) ) + ( m*cos(theta(1,a-1))*k2*theta(1,a-1) ) + (m*l*(theta(2,a-1))^2 *sin(theta(1,a-1))) + ( -k1*x(1,a-1) ));
        theta(3,a) = (-x(3,a)*cos(theta(1,a-1))/l) + ( sin(theta(1,a-1))*g/l) - (k2*theta(1,a-1)/l);
    
        x(2,a) = x(2,a-1) + x(3,a)*dt;
        x(1,a) = (x(1,a-1)+x(2,a-1)*dt + x(3,a)*(dt^2)*0.5);
    
        theta(2,a) = theta(2,a-1) + theta(3,a)*dt;
        theta(1,a) = (theta(1,a-1)+theta(2,a-1)*dt + theta(3,a)*(dt^2)*0.5);
    
        Time(a)=Time(a-1)+dt;
        if ( abs(theta(1,a)) > 1.5708  )
            boolean = 0;
        end
    end

    if (boolean == 0)
        boolean = 1;
    else
        break
    end
end

disp(k2);

figure(3)
 plot(Time(1,:) ,theta(1,:));
 title("theta againt time stage 3");
 xlabel("Time(s)");
 ylabel("Theta ");