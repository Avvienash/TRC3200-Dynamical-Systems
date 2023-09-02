%% Week 4: LAB 2
% Name: Avvienash A/L Jaganathan
% ID: 32281013
% Date: 25/3/2022
clear all; close all; clc;
%% Stage 1
% Question 1

a = pi/6; % alpha = 30 dg (convert to rad)
g = 9.81; % gravatational constant

v_th = 0.5; % constant speed of theta
% const speed = no acc

a_p = g*cos(a)*sin(a); % calculate acc of p

a_z = g; % acc along z

T = 15;
tsteps = 0.1;

mat=zeros(T/tsteps+1,5);
%initialise
mat(1,2) = 1000; % p_0 =1000
mat(1,3) = 1000 /(tan(a)); 
for loop = 1:(T/tsteps)
    
    % find for theta
    mat(loop+1,1) = rem((mat(loop,1) + v_th*tsteps),(2*pi));
    
    % find for z
    
    mat(loop+1,5)= mat(loop,5) - a_z*tsteps;
    mat(loop+1,3)= mat(loop,3) + mat(loop,5)*tsteps;
    
    % find for p
    mat(loop+1,4)= mat(loop,4) -a_p*tsteps;
    mat(loop+1,2)= mat(loop,2) +mat(loop,4)*tsteps;
end

figure(1)
subplot(2,1,1);
plot(mat(:,2),mat(:,4))
xlabel("p")
ylabel("dp")
title("dp against p")

subplot(2,1,2); 
plot(mat(:,3),mat(:,5))
title("dz against z")
xlabel("z")
ylabel("dz")

figure(2)
x = mat(:,2).*sin(mat(:,1));
y = mat(:,2).*cos(mat(:,1));
plot(x,y)
title("x agaisnt y")
xlabel("x")
ylabel("y")

%% Stage 2
% step 1 and 2
clear all; close all; clc;

a = pi/6; % alpha = 30 dg (convert to rad)
g = 9.81; % gravatational constant

z = 1000/tan(a);
t = sqrt(2*z/g)

b = pi/3; 

z2 = 1000/tan(b);
t2 = sqrt(2*z2/g)
%% Stage 2 Step 3
clear all; close all; clc;

a = pi/3; % alpha = 30 dg (convert to rad)
b = pi/6
% double alpha
g = 9.81; % gravatational constant

v_th = 0.5; % constant speed of theta
% const speed = no acc

a_p = g*cos(a)*sin(a); % calculate acc of p
a_2p = g*cos(b)*sin(b);

a_z = g; % acc along z

T = 11;
T2= 19;
tsteps = 0.1;

mat=zeros(T/tsteps+1,5);
%initialise
mat(1,2) = 1000; % p_0 =1000
mat(1,3) = 1000 /(tan(a)); 

mat2=zeros(T2/tsteps+1,5);
%initialise
mat2(1,2) = 1000; % p_0 =1000
mat2(1,3) = 1000 /(tan(b)); 
for loop = 1:(T/tsteps)
    
    % find for theta
    mat(loop+1,1) = rem((mat(loop,1) + v_th*tsteps),(2*pi));
    
    % find for z
    
    mat(loop+1,5)= mat(loop,5) - a_z*tsteps;
    mat(loop+1,3)= mat(loop,3) + mat(loop,5)*tsteps;
    
    % find for p
    mat(loop+1,4)= mat(loop,4) -a_p*tsteps;
    mat(loop+1,2)= mat(loop,2) +mat(loop,4)*tsteps;

end

for loop = 1:(T2/tsteps)
        
    %
    mat2(loop+1,1) = rem((mat2(loop,1) + v_th*tsteps),(2*pi));
    
    % find for z
    
    mat2(loop+1,5)= mat2(loop,5) - a_z*tsteps;
    mat2(loop+1,3)= mat2(loop,3) + mat2(loop,5)*tsteps;
    
    % find for p
    mat2(loop+1,4)= mat2(loop,4) -a_2p*tsteps;
    mat2(loop+1,2)= mat2(loop,2) +mat2(loop,4)*tsteps;
end

plot(mat(:,3),mat(:,2),mat2(:,3),mat2(:,2))
xlabel("z")
ylabel("p")
legend('2a','a','Location','southeast')
