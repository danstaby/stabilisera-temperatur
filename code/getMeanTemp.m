function coef = getMeanTemp(display)
%coef = getMeanTemp(display)
%
%Calculates the least square interpolation of the average
%temperature for Gothenburg the last twenty years.
%T = a + b*cos(omega*t) + c*sin(omega*t)
%with omega so it equals 2*pi for t = the number of days in a year.


if(nargin == 0)
  display = 1;
end

%Temperature data for the months
T = [0.7, 0.4, 2.6, 7.3, 11.9, 15.0, 17.8, 17.3, 13.4, 8.6, 4.6,1.4]';

d = (365-31)*[0:11]/11+15; %Create 


%Generate the system of equations
A = [ones(12,1), cos(2*pi*d/365)', sin(2*pi*d/365)'];

coef = A\T; %Solve least squares

dnew = 0:365;

if(display == 1) %Display data
  figure(1)
  plot(d,T, '*')
  hold on
  plot(dnew, coef(1)+coef(2)*cos(2*pi*dnew/365)+coef(3)*sin(2*pi*dnew/365), ...
       'r')
  hold off

  xlabel('Dag')
  ylabel('Medeltemperatur')
  xlim([0, 365])
  legend('Data', 'Interpolering')
end