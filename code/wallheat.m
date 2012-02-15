function ret = wallheat()
%wallheat()
%
%This script calculates the steady state temperatures and heat flows
%in a wall and plots the temperature dependance of the heat flow


BayArea = 47;
NorthArea = 290;
OtherArea = 151+61;

%BayWall is unknown. The previous report was not clear enough on
%this point to make us able to calculate this flow. 
	   
NorthWall = [0.02, 1; %Lengths and thermal conductivities have been taken
	    0.1, 0.037; %from the previous report.
	    0.02, 1;
	    0.5, 0.6;
	    0.01, 1];

OtherWall = [0.5, 0.6;
	    0.01, 1];

BoundaryTemperature = [20, 0];

HeatLoss = [];
T = [-30:0.1:30];
for currentTemp =-30:0.1:30
  
  BoundaryTemperature(2) = currentTemp;

  [HeatNorth, Temperature] = heatrod(BoundaryTemperature, NorthWall);
  [HeatOther, Temperature] = heatrod(BoundaryTemperature, OtherWall);
  if(currentTemp == -30)
    HeatNorth*NorthArea
  elseif(currentTemp == -29)
    HeatNorth*NorthArea
    NorthArea
  end
  %Add the walls adjacant to the bay windows.
  lossNorth = NorthArea*HeatNorth;
  lossOther = OtherArea*HeatOther;
  lossTotal = lossNorth + lossOther;
  HeatLoss = [HeatLoss; lossNorth, lossOther, lossTotal ];
end

plot(T, HeatLoss(:,1), 'b');
hold on
plot(T, HeatLoss(:,2), 'r');
plot(T, HeatLoss(:,3), 'k');
legend('North Wall', 'Other Walls', 'Total');
xlabel('Temperature outside in celsius')
ylabel('Heat loss (W)')

hold off