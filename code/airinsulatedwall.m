function ret = airinsulatedwall(StepLength, MaxLength)
%airinsulatedwall(StepLength, MaxLength)
%

TemperatureInside = 20;
TemperatureOutside = -20;
	   
Wall = [0.5, 0.6;
	0.01, 1];

LayerCount = size(Wall);
LayerCount = LayerCount(1);

BoundaryTemperature = [TemperatureInside, TemperatureOutside];

[WallHeat, Temperature] = heatrod(BoundaryTemperature, Wall);

InsulatedWall = [Wall; 0, 0.024];

Length = [0:StepLength:MaxLength];

AirHeat = [WallHeat];
Temp = [TemperatureOutside];
for L=StepLength:StepLength:MaxLength
  InsulatedWall(LayerCount + 1, 1) = L;
  [Heat, Temperature] = heatrod(BoundaryTemperature, ...
				InsulatedWall);
  AirHeat = [AirHeat Heat];
  Temp = [Temp Temperature(LayerCount)];
end

InsulatedWall
figure(1)
plot(Length, WallHeat,'k');
hold on
plot(Length, AirHeat);
xlabel('Air width (m)')
ylabel('Heat flow (W/m^2)')
hold off
figure(2)
plot(Length, 100*(1-AirHeat/WallHeat))
xlabel('Air width (m)')
ylabel('Maximum loss due to convection (%)')

figure(3)
plot(Length, Temp)