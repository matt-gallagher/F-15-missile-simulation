function [xdot] = dx(t,x)

g = [0; 0; -9.80665]; %gravitational accleration vector
d = .07; %diameter of rocket
r = d/2; %radius of rocket
s = pi * r^2; %area of rocket
pos_vec = x(1:3);
velo_vec = x(4:6);
[rho,a] = atmosphere_model(pos_vec);
M = norm(velo_vec)/a;

%% Table values

M_boost = [0.1 0.3019 0.4655 0.6295 0.8082 0.9984 1.1866 1.3931 1.6186 1.8782 2.0964 2.1066 2.25]'; % Mach number while boosting
drag_boost = [0.55 0.5502 0.5509 0.5505 0.5684 0.6841 0.7479 0.7555 0.7367 0.6925 0.6624 0.6611 0.64]'; % Drag coefficient while boosting
M_coast = [0.1 0.7841 0.8258 0.8616 0.9024 0.9500 1.0027 1.0956 1.2100 1.3371 1.4028 1.6031 1.8107 2.024 2.1066 2.25]'; % Mach number while coasting
drag_coast = [0.7 0.7031 0.7357 0.7711 0.8122 0.8791 0.9615 0.9988 1.0106 1.0053 0.9983 0.9647 0.9115 0.8543 0.834 0.8]'; % Drag coefficient while coasting
time_values = [0 .078 .1793 .2793 .3793 .4821 .5821 .6821 .7821 .8821 .9821 1.0821 1.112 1.113 100]'; %Timestamps after launch
thrust_values = [6920.5 6920.5 6153.2 5981.4 5858.1 5835.1 6132.4 6602.4 7072.5 7302 8760 2963.2 814.9 0 0]'; % Thrust values over time
mass_values = [10.4319 10.1886 9.8901 9.6176 9.3507 9.0666 8.7719 8.4845 8.1759 7.8474 7.4813 7.1776 7.1545 7.1368 7.1368]'; % Rocket mass over time

%% Transition time at 1.113 sec

if t > 1.113
    Cd = interp1(M_coast, drag_coast, M, 'linear', 'extrap'); %Cd at time t during coast
else
    Cd = interp1(M_boost, drag_boost, M, 'linear', 'extrap'); %Cd at time t duting boost
end

mass = interp1(time_values, mass_values, t, 'linear', 'extrap'); %mass at time t
thrust = interp1(time_values, thrust_values, t, 'linear', 'extrap'); %thrust at time t

%% Acceleration

thrust_accel = (thrust*velo_vec/norm(velo_vec))/mass; %acceleration due to thrust
drag_accel = (-.5*rho*norm(velo_vec)*Cd*s.*velo_vec)/mass; %acceleration due to drag
acceleration = thrust_accel + drag_accel + g; %total acceleration
xdot = [velo_vec; acceleration];
end