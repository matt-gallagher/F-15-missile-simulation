clc
clear
format long %clears command window and formats output
%% Initial Conditions

psi = 0; % azimuth angle
theta = -30; %elevation angle
Mo = 0.5; %initial mach number
N = 5; %numbr of stochastic runs

sigmapsi = 1; %standard deviation psi
sigmatheta = 1; %standard deviation theta
sigmaM = .02; %standard deviation Mach

to = 0; %initial time
tf = 100; % final time
dt = .01; %time step
steps = tf/dt; %number of steps
Ztrue = 100; %impact altitude
pos = [0; 0; 1500]; % initial position
Mach(1) = Mo; %save first value of Mach number
[rho, a] = atmosphere_model(pos); %atmposhere model
speed = a * Mo; %initial speed
Vo = speed*[cosd(theta)*sind(psi); cosd(theta)*cosd(psi); sind(theta)]; %velocity in vector form
X_state = [pos; Vo]; %state vector

[xdot] = dx(to, X_state); %xdot of initial conditions

velo(1) = norm(X_state(4:6)); %save first speed value
output(1,:) = [to, X_state', xdot(4:6)']; %output initial time, position, velocity, and acceleration

%% Deterministic Run

t = to; %start at initial time
for ind = 2:steps
    [time, x] = ode45(@dx, [t t+dt], X_state); %vlaues of X_state at t+dt
    t = time(end); %last time value
    X_state = x(end, :); %state at last time value
    [xdot] = dx(t,X_state'); %velocity and acceleration from xdot
    output(ind, :) = [t, X_state, xdot(4:6)']; %output time, position, velocity, acceleration
    velocity(ind) = norm(X_state(4:6)); %output speed
    Mach(ind) = velocity(ind)./a; %compute Mach number
    
    %termination condition
    if X_state(3) < Ztrue %termination condition
       t_final = (Ztrue - output(end, 4))/output(end, 7); %time before impact
       t_impact = output(end, 1) + t_final; %time of impact
       Impact = (output(end, 2:4) + output(end, 5:7) * t_final); %point of impact
       Data = [t_impact, Impact]; %vector with time and position of impact
       break;
    end
    
end

%% Outputs

Xtrue = Data(2);
Ytrue = Data(3);

disp('Scenario A')
fprintf('Impact time = %5.2f seconds \n', Data(1)) % displays time
fprintf('Impact x position = %5.2f meters \n', Xtrue) % displays x impact location
fprintf('Impact y position = %5.2f meters  \n', Ytrue) % displays y impact location
fprintf('Impact z position = %5.2f meters  \n', Data(4))% displays z impact location
 
% speed vs. time 
figure(1)
plot(output(:,1), velocity, 'linewidth',2)
grid on;
title('Speed vs. time')
xlabel('time')
ylabel('Speed')
hold on

% mach vs. time
figure(2)
plot(output(:,1), Mach, 'linewidth',2)
grid on;
title('Mach number vs. time')
xlabel('time')
ylabel('Mach number') 

%% Monte Carlo

for i = 1:N
    Mach_rand = sigmaM*randn; %random Mach values
    theta_rand = sigmatheta*randn; %random theta values
    psi_rand = sigmapsi*randn; %random psi values
    theta_st = theta + theta_rand;
    psi_st = psi + psi_rand;
    Mach_st = Mo + Mach_rand;
[rho, a] = atmosphere_model(pos); %atmposhere model
speed_st = a * Mach_st; %initial speed
Vo_st = speed_st*[cosd(theta_st)*sind(psi_st); cosd(theta_st)*cosd(psi_st); sind(theta_st)]; %velocity in vector form
X_state_st = [pos; Vo_st]; %state vector

[xdot] = dx(to, X_state_st); %xdot of initial conditions

velo(1) = norm(X_state_st(4:6)); %save first speed value
output_st(1,:) = [to, X_state_st', xdot(4:6)']; %output initial time, position, velocity, and acceleration

t = to; %start at initial time
for ind = 2:steps
    [time, x] = ode45(@dx, [t t+dt], X_state_st); %vlaues of X_state at t+dt
    t = time(end); %last time value
    X_state_st = x(end, :); %state at last time value
    [xdot] = dx(t,X_state_st'); %velocity and acceleration from xdot
    output_st(ind, :) = [t, X_state_st, xdot(4:6)']; %output time, position, velocity, acceleration
    velocity_st(ind) = norm(X_state_st(4:6)); %output speed
    Mach(ind) = velocity_st(ind)./a; %compute Mach number
    
    %termination condition
    if X_state_st(3) < Ztrue %termination condition
       t_final = (Ztrue - output_st(end, 4))/output_st(end, 7); %time before impact
       t_impact = output_st(end, 1) + t_final; %time of impact
       Impact = (output_st(end, 2:4) + output_st(end, 5:7) * t_final); %point of impact
       Data = [t_impact, Impact]; %vector with time and position of impact
       dX(i) = Data(2) - Xtrue
       dY(i) = Data(3) - Ytrue
       break;
    end
    
end

end
mu_x = mean(dX);
mu_y = mean(dY);
std_x = std(dX);
std_y = std(dY);
disp(' ')
disp(' ')
disp('Stochastic Simulation')
fprintf('Standard Deviation x = %5.2f \n', std_x)
fprintf('Standard Deviation y = %5.2f \n', std_y)
fprintf('Mean Deviation x = %5.2f \n', mu_x)
fprintf('Mean Deviation y = %5.2f \n', mu_y)


