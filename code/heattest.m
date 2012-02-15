function ret = heattest()
%heattest()
%
%This script is a test file for the heatrod script. The values
%used in this script is chosen to compare the scripts output
%to manually calculated values

Layers = [0.1, 100; %This is the layer vector. The first column 
	 0.15, 15; %states the width of the layer and the second 
	  0.4, 80]; %column states the thermal conductivity.

BoundaryTemperature = [400, 100];

[Heat, Temperature] = heatrod(BoundaryTemperature, Layers)