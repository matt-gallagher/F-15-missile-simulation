function [rho,a] = atmosphere_model(pos)
%% Constants

To = 288.15; %temp at MSL
Po = 101325; %pressure at MSL
L = -0.0065; %lapse rate
zo = 0; %MSL altitude
g = 9.80665; %gravitational acceleration
R = 287.05; %J/kg-K
gam = 1.4; %specific heat ratio

%% Model

z = pos(3); %z is third term in position vector
T = To + L * (z - zo); %Temp at current altitude
gLR = -(g/(L * R)); %Pressure equation exponent
P = Po * (T/To)^gLR; %Pressure at curernt altitude
rho = P/(R * T); %Density at current altitude
a = sqrt(gam * R * T); %Acoustic speed
end
