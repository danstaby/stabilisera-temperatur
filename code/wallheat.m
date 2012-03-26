function ret = wallheat()
%wallheat()
%
%This script calculates the steady state temperatures and heat flows
%in a wall and plots the temperature dependance of the heat flow


BayArea = 47;
NorthArea = 290;
OtherArea = 151+61;

	   
NorthWall = [0.02, 1; %Lengths and thermal conductivities have been taken
	    0.1, 0.037; %from the previous report.
	    0.02, 1;
	    0.5, 0.6;
	    0.01, 1];

BayWall = [0.024, 0.25;
	   0.05, 0.037;
	   0.016, 0.25;
	   0.025, 0.025;
	   0.016, 0.5;
	   0.002, 401];

OtherWall = [0.5, 0.6;
	    0.01, 1];

BoundaryTemperature = [20, 0];

HeatLoss = [];
T = [-30:0.1:30];

for currentTemp =-30:0.1:30
  
  BoundaryTemperature(2) = currentTemp;

  [HeatNorth, Temperature] = heatrod(BoundaryTemperature, NorthWall);
  [HeatOther, Temperature] = heatrod(BoundaryTemperature, OtherWall);
  [HeatBay, Temperature] = heatrod(BoundaryTemperature, BayWall);
  
  lossNorth = NorthArea*HeatNorth;
  lossOther = OtherArea*HeatOther;
  lossBay = BayArea*HeatBay;
  lossTotal = lossNorth + lossOther + lossBay;
  HeatLoss = [HeatLoss; lossNorth, lossOther, lossBay,  lossTotal ];
end

plot(T, HeatLoss(:,4), 'k');
xlabel('Temperature outside in celsius')
ylabel('Heat loss (W)')

figure(2)
plot(T,HeatLoss(:,2)/OtherArea,'k')


disp(['Loss from the wall adjacant to the bay windows: '...
     num2str(100*HeatLoss(1,3)/HeatLoss(1,4)) '%']);
disp(['Loss from the north wall: ' ...
      num2str(100*HeatLoss(1,1)/HeatLoss(1,4)) '%']);
disp(['Loss from the other walls: ' ...
     num2str(100*HeatLoss(1,2)/HeatLoss(1,4)) '%']);