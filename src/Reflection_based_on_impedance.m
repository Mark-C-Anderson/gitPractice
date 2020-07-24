% Purpose: to understand the amount of energy that is reflected when sound
% waves hit a barrier between two acoustic impedances. This is to help with
% the results presented in the Characterization of the Effects of Ground
% Boards on Acoustic Signals research that I did as an intern at NASA
% Langley Research Center in Summer 2019.

% including some nice stuff. This is no longer needed in the new file structure
addpath('/Users/markanderson/Box/Acoustics MCA/MATLAB Code/Simple Analysis Functions/')
addpath('/Users/markanderson/Box/Acoustics MCA/MATLAB Code/Plots/')
setPlotStyle('Group Standard')

% List the characteristic impedances
% air = 1.225 * 343 = 420
% sand = 1600 * 1626 = 2.60e6 (google)
% wood = 680 * 4000 = 2.76e6 (google and Engineering Toolbox)
% plastic = 1200 * 1000 = 1.2e6 (google)
Z1 = 1.2e6;
Z2 = 2.6e6;

% Give incident angle range in degrees
incidentAngle = 0:0.2:90;
incidentAngle = incidentAngle*pi/180; % Converting to radians for calculations

% Solve for the speeds of sound
rho_0 = 101325; % atmospheric pressure
c1 = Z1/rho_0; % eqn. C-60 (Blackstock)
c2 = Z2/rho_0; % eqn. C-60 (Blackstock)

% Solve for the transmission angle
transmissionAngle = asin(c2/c1.*sin(incidentAngle)); % Snell's Law

% Solve for the reflectance coefficient
R = (Z2.*cos(incidentAngle) - Z1.*cos(transmissionAngle)) ./ ...
    (Z2.*cos(incidentAngle) + Z1.*cos(transmissionAngle));

% Solve for the transmission coefficient
T = 2*Z2.*cos(incidentAngle) ./ (Z2.*cos(incidentAngle) + ...
    Z1.*cos(transmissionAngle));

% Plotting the transmission angle
figure()
plot(incidentAngle*180/pi,transmissionAngle*180/pi)
title('Transmission Angle')
xlabel('Incident Angle (Degrees (re: normal))')
ylabel('Tranmission Angle (Degrees (re: normal))')
grid on

% Plotting the reflectance coefficient
figure()
plot(incidentAngle*180/pi,R)
hold on
plot(incidentAngle*180/pi,T)
plot(incidentAngle*180/pi,T-R)
title('Reflection and Transmission Coefficients')
xlabel('Incident Angle (Degrees (re: normal))')
ylabel('Coefficient Value')
grid on
legend('Reflectance','Transmission','T-R = 1')