function ret = airinsulatedwall(AirWidth)
%airinsulatedwall(AirWidth)
%
%This script calculates the steady state temperatures in an air
%insulated wall for the WallSimulation application

if(nargin == 0)
  AirWidth = 1.5;
end


	   
Wall = [0.5, 0.6;
	0.01, 1;
	0.024, AirWidth];

BoundaryTemperature = [20, 0];

WallTemp = [];
T = [-30:0.1:30];
for currentTemp =-30:0.1:30
  
  BoundaryTemperature(2) = currentTemp;

  
  [Heat, Temperature] = heatrod(BoundaryTemperature, Wall);
  WallTemp = [WallTemp Temperature(max(size(Temperature)))];
end

plot(T, WallTemp, 'k');

xlabel('Temperature outside in celsius')
ylabel('Wall temperature')

