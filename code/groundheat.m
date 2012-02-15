function ret = groundheat()
%This scripts purpose is to calculate the heatflow through the
%foundation of the building.

Layers = [0.1, 100; %This is the layer vector. The first column 
	 0.15, 15; %states the width of the layer and the second 
	  0.4, 80]; %column states the thermal conductivity.

BoundaryTemperature = [400, 100];

[Heat, Temperature] = heatrod(BoundaryTemperature, Layers)