function ret = groundheat()
%This scripts purpose is to calculate the heatflow through the
%foundation of the building.

Layers = [0.7, 1; %0.25m concrete. (Width? Type? Heat conductivity?) 
	 2, 0.024; %2m air (Width?)
	  0.7, 1; %0.25m concrete (Type? Heat conductivity?)
	  0.3, 0.024; %Air
	  1.5, 2.2]; %Granite

Area = 285;

BoundaryTemperature = [20, 6];

[Heat, Temperature] = heatrod(BoundaryTemperature, Layers);
Temperature
HeatLoss = Heat*Area