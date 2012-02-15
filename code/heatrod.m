function [Heat, Temperature] = heatrod(BoundaryTemperature,Layers)
%This script calculates the steady state temperatures and
%heat flows in a rod.
%
%[Heat, Temperature] = heatrod(BoundaryTemperature, Layers)
%
%Inputs:
%BoundaryTemperature - Vector with boundary temperatures in centigrades
%Layers - Matrix with length of rod elements in first column and
%thermal conductivity in the second column
%
%Outputs:
%Heat - The heat flow per surface unit (W/m^2)
%Temperature - Vector with the internal temperatures in centigrades

kelvin = 273.15; %0 celsius in kelvin

LayerCount = size(Layers);
LayerCount = LayerCount(1);

BoundaryTemperature = kelvin + BoundaryTemperature;

fourierHeat = [-1, 1; 1, -1]; %Seed for stiffness matrix

%Create stiffness matrix

stiffness = zeros(LayerCount+1);
for n=1:LayerCount
  stiffness(n:(n+1), n:(n+1)) = stiffness(n:(n+1),n:(n+1)) - ...
    Layers(n,2)*fourierHeat/Layers(n,1);
end

%Multiply boundary values to the appropriate columns in the
%stiffness matrix

stiffness(:,1) = BoundaryTemperature(1)*stiffness(:,1);
stiffness(:,LayerCount+1) = BoundaryTemperature(2)*stiffness(:, ...
						  LayerCount+1);

%Move the unknown boundary heat flows to the x-side of the equation
%and move the first and last column to the y side of the equation

y = -stiffness(:,1)-stiffness(:,LayerCount+1);

%Create new stiffness matrix 

A = zeros(LayerCount+1);
A(1,1) = -1;
A(LayerCount+1, LayerCount+1) = -1;

A(1:LayerCount+1, 2:LayerCount) = ...
    stiffness(1:LayerCount+1, 2:LayerCount);

%Solve equation

x = A\y;
Heat = x(1);
Temperature = x(2:LayerCount) - kelvin;
