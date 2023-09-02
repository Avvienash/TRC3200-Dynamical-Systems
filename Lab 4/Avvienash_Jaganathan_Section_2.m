%% Week 6: LAB 4
% Name: Avvienash A/L Jaganathan
% ID: 32281013
% Date: 29/4/2022
clear all; close all; clc;

m = 1;
l = 3;
g = 9.81;
k = 2;
R = 0.5;
dt = 2e-4;
n =10/(2e-4);
t = 0:dt:dt*(n-1);
r = 5:20;
timetaken = zeros(1,16);

for b = 1:16
    q = zeros(3,n);
    q(2,1) = 8;
    speed = ones(1,10);
    rf = r(b);

    for a = 2:n
    
        % use eq motion to find ac
        term1 = m*l*l;
        term2 = m*g*l*sin(q(1,a-1));
        term3 = 4*k*(R^2)*sin(q(1,a-1))*cos(q(1,a-1));
        term4 = rf*q(2,a-1);
        q(3,a)= -(term2+term3+term4)/term1;
        
        % use a to calc rest
        q(2,a) = q(2,a-1) + q(3,a)*dt;
        q(1,a) = (q(1,a-1)+q(2,a-1)*dt + q(3,a)*(dt^2)*0.5);

        speed = circshift(speed,1);
        speed(1) = q(2,a);

        bool = (abs(speed(1,:)) <= (1e-4));
        if all(bool == 1)
            timetaken(1,b) = a;
        else
            continue
        end
        break
    end
 
end

Min = min(timetaken(timetaken>0));
Ind = find( timetaken == Min);

RF = r(Ind);
TIME = Min*dt;
fprintf("The time is shortest at Rf = %d and the time taken is %d \n",RF,TIME)

%%
fprintf("R double : \n")
clear all

m = 1;
l = 3;
g = 9.81;
k = 2;
R = 1;
dt = 2e-4;
n =10/(2e-4);
t = 0:dt:dt*(n-1);
r = 5:20;
timetaken = zeros(1,16);

for b = 1:16
    q = zeros(3,n);
    q(2,1) = 8;
    speed = ones(1,10);
    rf = r(b);

    for a = 2:n
    
        % use eq motion to find ac
        term1 = m*l*l;
        term2 = m*g*l*sin(q(1,a-1));
        term3 = 4*k*(R^2)*sin(q(1,a-1))*cos(q(1,a-1));
        term4 = rf*q(2,a-1);
        q(3,a)= -(term2+term3+term4)/term1;
        
        % use a to calc rest
        q(2,a) = q(2,a-1) + q(3,a)*dt;
        q(1,a) = (q(1,a-1)+q(2,a-1)*dt + q(3,a)*(dt^2)*0.5);

        speed = circshift(speed,1);
        speed(1) = q(2,a);

        bool = (abs(speed(1,:)) <= (1e-4));
        if all(bool == 1)
            timetaken(1,b) = a;
        else
            continue
        end
        break
    end
 
end

Min = min(timetaken(timetaken>0));
Ind = find( timetaken == Min);

RF = r(Ind);
TIME = Min*dt;
fprintf("The time is shortest at Rf = %d and the time taken is %d \n",RF,TIME)

%%
fprintf("k double : \n")
clear all

m = 1;
l = 3;
g = 9.81;
k = 4;
R = 0.5;
dt = 2e-4;
n =10/(2e-4);
t = 0:dt:dt*(n-1);
r = 5:20;
timetaken = zeros(1,16);

for b = 1:16
    q = zeros(3,n);
    q(2,1) = 8;
    speed = ones(1,10);
    rf = r(b);

    for a = 2:n
    
        % use eq motion to find ac
        term1 = m*l*l;
        term2 = m*g*l*sin(q(1,a-1));
        term3 = 4*k*(R^2)*sin(q(1,a-1))*cos(q(1,a-1));
        term4 = rf*q(2,a-1);
        q(3,a)= -(term2+term3+term4)/term1;
        
        % use a to calc rest
        q(2,a) = q(2,a-1) + q(3,a)*dt;
        q(1,a) = (q(1,a-1)+q(2,a-1)*dt + q(3,a)*(dt^2)*0.5);

        speed = circshift(speed,1);
        speed(1) = q(2,a);

        bool = (abs(speed(1,:)) <= (1e-4));
        if all(bool == 1)
            timetaken(1,b) = a;
        else
            continue
        end
        break
    end
 
end

Min = min(timetaken(timetaken>0));
Ind = find( timetaken == Min);

RF = r(Ind);
TIME = Min*dt;
fprintf("The time is shortest at Rf = %d and the time taken is %d \n",RF,TIME)

%%
fprintf("Double m : \n")
clear all

m = 2;
l = 3;
g = 9.81;
k = 2;
R = 0.5;
dt = 2e-4;
n =10/(2e-4);
t = 0:dt:dt*(n-1);
r = 5:20;
timetaken = zeros(1,16);

for b = 1:16
    q = zeros(3,n);
    q(2,1) = 8;
    speed = ones(1,10);
    rf = r(b);

    for a = 2:n
    
        % use eq motion to find ac
        term1 = m*l*l;
        term2 = m*g*l*sin(q(1,a-1));
        term3 = 4*k*(R^2)*sin(q(1,a-1))*cos(q(1,a-1));
        term4 = rf*q(2,a-1);
        q(3,a)= -(term2+term3+term4)/term1;
        
        % use a to calc rest
        q(2,a) = q(2,a-1) + q(3,a)*dt;
        q(1,a) = (q(1,a-1)+q(2,a-1)*dt + q(3,a)*(dt^2)*0.5);

        speed = circshift(speed,1);
        speed(1) = q(2,a);

        bool = (abs(speed(1,:)) <= (1e-4));
        if all(bool == 1)
            timetaken(1,b) = a;
        else
            continue
        end
        break
    end
 
end

Min = min(timetaken(timetaken>0));
Ind = find( timetaken == Min);

RF = r(Ind);
TIME = Min*dt;
fprintf("The time is shortest at Rf = %d and the time taken is %d \n",RF,TIME)