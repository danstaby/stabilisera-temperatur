function ret = walltransientfemdecember(Time, Nodes, Insulated, UseSun, StaticTemp)
% walltransientfem(Time, Nodes, Insulated, UseSun, StaticTemp)
%
% This Script calculates the transient temperatures in a wall with
% a specified amount of nodes. Insulated = 0 means uninsulated wall
% and Insulated = 1 mean that the wall is insulated.
%StaticTemp = -999 means that the dynamic temperature is used instead.
%UseSun = 0 means no sun, UseSun=1 means full sun, UseSun = 2 means cloudy

global bUseSun
bUseSun = UseSun;

%Temperature inside

PrepareInterpolation('sunintensity_december.txt');
PrepareTempInterpolation([6, -11;16,-5]);

%t = [0:60:24*3600];
%temp = zeros(1,size(t,2));
%sun = zeros(1,size(t,2));
%for n = 1:size(t,2)
  
%  temp(n) = Tout(t(n));
%  sun(n) = Qsun(t(n));
%end

%figure(1)
%hold off
%plot(t/3600, temp)
%figure(2)
%plot(t/3600, sun)



Tin = 20;

if(nargin < 5)
  StaticTemp = -999;
end
if(nargin < 4)
  UseSun = 1;
end
if(nargin < 3)
  Insulated = 0;
end


if( Insulated == 0)
  Material = [0.5, 5.2e-7, 0.6, 1153846.15];
  % Material = [Width, Thermal Diffusivity (depricated), Thermal
  % Conductivity, Volumetric heat capacaty (density times specific heat)]
  

else
  Material = [0.1, 0.037/(1400*840), 0.037, 840*70;
	      0.5, 5.2e-7, 0.6, 1153846.15];

  %110 kg/m^3
  %840 J/K kg
  %0.037 W/mK
end


kelvin = 273.15; %0 centigrades in kelvin


%Calculate the total width of the wall (the sum of the
%constitutant parts)

Width = sum(Material(:,1));

%Initialize some variables
NrNodes = floor(Nodes*Material(:,1)/Width)+1;

NrNodes(NrNodes < 2)  = 2; %There needs to be at least one interior
                           %point in each piece of wall.

MaterialCount = size(Material,1); %Number of different materials in
                                  %the wall.
NodeCount = sum(NrNodes)+1; %Count of all nodes
Pts = zeros(NodeCount,1); %Vector that contains all the discrete nodes
L = zeros(NodeCount,1); %Length of elements
uInitial = sparse(NodeCount,1); %Initial Guess
rho = zeros(NodeCount-1,1); %Volumetric heat capacity
kappa = zeros(NodeCount-1,1); %Thermal conductivity
PrevEnd = 0;
PrevL = 0;

%Prepare the discretization of the spatial dimension
for n = 1:size(Material,1)
  Pts(PrevEnd+(1:NrNodes(n))) = PrevL + Material(n,1)*(0: ...
						  (NrNodes(n)-1))'/(NrNodes(n));
  rho(PrevEnd+(1:NrNodes(n))) = Material(n,4);
  kappa(PrevEnd+(1:NrNodes(n))) = Material(n,3);
  L(PrevEnd + (1:NrNodes(n))) = Material(n,1)/(NrNodes(n)); %length
                                                            %of the
                                                            %wall pieces
  PrevL = PrevL + Material(n,1); %Ad hoc iteration variables
  PrevEnd = PrevEnd + NrNodes(n);
end

Pts(sum(NrNodes)+1) = PrevL;


%Set initial value to the steady state solution
if(size(Material,1) > 1)
  [~, TempInit] = heatrod([-8 Tin], Material(:, [1 3]));
end
PrevEnd = 0;
Tlast = 0;

for n = 1:(MaterialCount-1)
  
  uInitial(PrevEnd+(1:NrNodes(n))) = (Pts(PrevEnd+(1:NrNodes(n)))-sum(L(1:(PrevEnd))))* ...
                                    (TempInit(n)-Tlast)/Material(n,1)+Tlast;
  Tlast = TempInit(n);
  PrevEnd = PrevEnd + NrNodes(n);
end

uInitial(PrevEnd+(1:(NrNodes(MaterialCount)))) = (Pts(PrevEnd+(1:(NrNodes(MaterialCount))))-...
					   sum(L(1:(PrevEnd))))*(Tin-Tlast)/Material(MaterialCount,1)...
	                                   +Tlast;
uInitial(NodeCount) = Tin;
uInitial = uInitial + kelvin;
uLast = uInitial;


%Initialization of the matrices needed to do the FEM calculations
A = sparse(NodeCount, NodeCount); %Stiffness matrix
Q = sparse(NodeCount, NodeCount); %Neumann contribution to
                                  %stiffness matrix
M = sparse(NodeCount, NodeCount); %Mass matrix
u = sparse(NodeCount,1); %Temperature vector
f = sparse(NodeCount,1); %Load vector
b = sparse(NodeCount,1); %Intermediary load vector (used to move known
                     %terms to the RHS)
g = sparse(NodeCount,1); %Neumann contribution to the load vector

Vec = zeros(NodeCount,NodeCount); %Eigen vectors
lambda = zeros(Nodes,1); %Eigennumbers of the problem (inverse of
                         %mass matrix times the stiffness matrix)
Energy = zeros(floor((Time(3)-Time(1))/Time(2))+1,1);
EnergyIn = zeros(floor((Time(3)-Time(1))/Time(2))+1,1);

dirbc = [NodeCount; %Column vectors with the dirichlet condition(s)
	Tin+kelvin];

%Assemble stiffness matrix

for n=1:(NodeCount-1)
  A([n n+1], [n n+1]) = A([n n+1], [n n+1]) + kappa(n)*[1/L(n),-1/L(n); -1/L(n), 1/L(n)];
end

%Assemble mass matrix

for n = 1:(NodeCount-1)
  M([n n+1], [n n+1]) = M([n n+1], [n n+1]) + rho(n)*L(n)*[1/3, 1/6; 1/6, 1/3];
end


Free = 1:(NodeCount-1); %Vector of the free nodes (i.e. the nodes
                        %without dirichlet conditions)

Const = sparse(NodeCount,1); %Intermidiary solution vector

Minv = sparse(NodeCount,NodeCount); %Inverse version of mass matrix
MinvA = sparse(NodeCount,NodeCount);%to speed up the computation

ab = sparse(NodeCount,1); %Steady-state solution

figure(1)
%plot(Pts, uInitial-kelvin)
hold on 

%Calculate the eigenvalues and the eigenvectors for the problem.
%These will be used to solve the system of linked ODE's.



[VecTemp lambdaTemp] = eig(full(inv(M(Free,Free))*A(Free,Free)));
Vec(Free,Free) = VecTemp; 
lambda(Free) = diag(lambdaTemp);


Minv(Free,Free) = inv(M(Free,Free));
MinvA(Free,Free) = inv(M(Free,Free))*A(Free,Free);


n = 0;
sigma = 5.670e-8;
Tguess = Tout(0) + kelvin;
Tolerance = 1e-3;
for t = Time(1):Time(2):Time(3)
  n = n+1;
  
  %Update boundary conditions
  To = StaticTemp;

  if(To == -999)
    To = Tout(t)+kelvin;
  end
  
  Qd = UseSun;
  Qs = 0;
  if(Qd ~= 0)
    [Qw, Qd] = Qsun(t);
  end
  
  Qtot = Qw + 0.2*Qd;
  h = 35;
  Rair = sigma*To^4*(1-0.261*exp(-7.77e-4*(273-To)^2));
  Raimb = sigma*To^4;
  Rtot = 0.4*Raimb + 0.6*Rair;
  

  Q(1,1) = -(-h); %T dependent neumann conditions
  g(1,1) = Qtot + h*To + Rtot - sigma*uLast(1)^4; %Constant neumann conditions

  %Move dirichlet conditions to the RHS
  u(:) = 0;
  u(dirbc(1,:)) = dirbc(2,:)';
  b = f - A*u + g;
  
  %Calculate eigenspace and solve the equations
  
  [VecTemp lambdaTemp] = eig(full(MinvA(Free,Free)+Minv(Free,Free)*Q(Free,Free)));
  Vec(Free,Free) = VecTemp;
  lambda(Free) = diag(lambdaTemp);

  ab(Free) = inv(A(Free,Free)+Q(Free,Free))*b(Free); %Ainv(Free,Free)*b(Free);
  Const(Free) = inv(Vec(Free,Free))*((uLast(Free)-ab(Free))); %./lambda(Free))
  u(Free) = Vec(Free,Free)*(exp(-lambda(Free)*Time(2)).*Const(Free)) + ab(Free);
  Tres = abs(Tguess-u(1));
  
  %fprintf(1, '\n');
  uLast(:) = u;
  
  %if(mod(t,1000) == 0)
  %uLast
  plot(Pts, uLast-kelvin, 'r');
  %disp('Yay plotted!')
  %end
  
  Energy(n) = kappa(NodeCount-1)*(uLast(NodeCount)-uLast(NodeCount- ...
						  1))/(L(NodeCount-1));
  EnergyIn(n) = kappa(1)*(uLast(1)*Q(1,1)+g(1,1));
end


hold off
figure(2)
plot((Time(1):Time(2):Time(3))/(24*3600), Energy)
hold on
plot((Time(1):Time(2):Time(3))/(24*3600)-1.5, tsmovavg(Energy,'s', floor(3*24* ...
						  3600/Time(2))+1,1), ...
     'r');

xlabel('Time (days)')
ylabel('Energy loss (W/m^2)')
legend('Momentary energy loss', 'Moving average (3 days)')
hold off


ret = [(Time(1):Time(2):Time(3))', Energy, EnergyIn];

function ret = PrepareInterpolation(fileName)
global SunResolution SplineInterpolation

%Load file
data = importdata(fileName);
time = data(:,2);
I = data(:,1);
SunResolution = mean(diff(time));
timesteps = size(time,1);

%Set up the system of equations
SplineInterpolation = zeros(2*(timesteps-1),1);
A = zeros(2*(timesteps-1), 2*(timesteps-1));
b = zeros(2*(timesteps-1),1);
for n = 1:(timesteps-1)
  A([2*n-1 2*n], [2*n-1 2*n]) = [ones(2,1), time([n n+1])];
  b([2*n-1 2*n]) = I([n n+1]);
end


%Solve for X 
SplineInterpolation = A\b;


function ret = PrepareTempInterpolation(points)
global BreakTimes TemperatureSpline

%Points is a matrix that contains the known
%times in the first column and the assosciated temperatures
%in the second column. The function should be periodic on 24h.

splineCount = size(points,1)+1;
A = zeros(splineCount*2, splineCount*2);
b = zeros(splineCount*2,1);
BreakTimes = points(:,1);

%Prepare equation system

for n=1:(splineCount-1)
  A(2*n-1, [2*n-1, 2*n]) = [1, points(n,1)];
  A(2*n, [2*(n+1)-1, 2*(n+1)]) = [1, points(n,1)];
  b([2*n-1 2*n]) = points(n,2);
end

A(2*splineCount-1,[1,2,2*splineCount-1, 2*splineCount]) = ...
    [1, 0,-1,-24];

A(2*splineCount,[2, 2*splineCount]) = [1, -1];

%Solve for X.
TemperatureSpline = A\b;

function [Qproj, Qdirect] = Qsun(Time)
global SunResolution SplineInterpolation bUseSun

t = mod(Time,24*3600);

if(((t/3600)-2) >= 0)
  [height, angle] = sunposition(11.979435,57.691522,2011,04,15,t/ ...
				(3600)-2);
else
  angle = 0;
  height = -1;
end

theta = mod(angle - (58+90),360);

projection = cosd(theta)*cosd(height);    

if(projection < 0 || height < 0)
  projection = 0;
end

spline = floor(t/(3600*SunResolution))+1; %Calculate which spline that
                                %should be used


I = [1, t/3600]*SplineInterpolation([2*spline-1 2*spline]);

%Get the solar intensity normal to the sun rays.
if(bUseSun == 1)
  Qproj = I*projection;
  Qdirect = I;
elseif(bUseSun == 2)
  Qproj = 0.2*I;
  Qdirect = 0.2*I;
end


function ret = Tout(Time)
global BreakTimes TemperatureSpline

t = mod(Time,24*3600)/3600;


foundSpline = 0;

for n = 1:(size(BreakTimes,1)-1)
  
  if((t > BreakTimes(n)) & (t < BreakTimes(n+1)))
    foundSpline = n+1;
  end
end

if(foundSpline == 0)
  
  if(t<= BreakTimes(1))
    foundSpline = 1;
  else
    foundSpline = size(BreakTimes,1)+1;
  end
end


ret = [1,t]*TemperatureSpline(2*foundSpline+[-1:0]);

%ret = 5*sin(pi*t/43200);
%T(t)=3/10*t+4.2, då 06.00<t<16.00
%T(t)=3/14*t+39/7 då t<06.00, t>16.00