function ret = walltransientfemstep(Time, Nodes, Insulated, Tstart,Tend)
% walltransientfemstep(Time, Nodes, Insulated, Tstart, Tend)
%
% This Script calculates the transient temperatures in a wall with
% a specified amount of nodes. Insulated = 0 means uninsulated wall
% and Insulated = 1 mean that the wall is insulated.
%Tstart is the start (steady state) temperature
%and Tend is the final temperature. These are the known
%temperatures of the outer boundary. The temperatuer inside the
%building is fixed to 20 centigrades



Tin = 20;


if( Insulated == 0)
  Material = [0.5, 5.2e-7, 0.6, 1153846.15];
  % Material = [Width, Thermal Diffusivity (depricated), Thermal
  % Conductivity, Volumetric heat capacaty (density times specific heat)]
  

elseif( Insulated == 1)
  Material = [0.1, 0.037/(1400*840), 0.037, 840*70;
	      0.5, 5.2e-7, 0.6, 1153846.15];

  %110 kg/m^3
  %840 J/K kg
  %0.037 W/mK
elseif(Insulated == 2)
  Material = [0.02, 0, 401, 390*8940;
	      0.016, 0, 0.5, 1.3e6;
	      0.016, 0, 0.25, 1090*800;
	      0.05, 0, 0.037, 840*70;
	      0.024, 0 ,0.25, 1090*800];

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
  [~, TempInit] = heatrod([Tstart Tin], Material(:, [1 3]));
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

  
%Update boundary conditions
To = Tend + kelvin;

h = 25;
Q(1,1) = -(-h); %T dependent neumann conditions
g(1,1) = (h*To); %Constant neumann conditions

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
n = 0;
figure(1)
hold on
for t = Time(1):Time(2):Time(3)
  n = n + 1;
  
  u(Free) = Vec(Free,Free)*(exp(-lambda(Free)*t).*Const(Free)) + ab(Free);
  
  
  plot(Pts, u-kelvin, 'r')
  Energy(n) = kappa(NodeCount-1)*(u(NodeCount)-u(NodeCount-1))/(L(NodeCount-1));
  EnergyIn(n) = kappa(1)*(u(1)*Q(1,1)+g(1,1));

end


hold off
figure(2)
plot((Time(1):Time(2):Time(3))/(24*3600), Energy)

xlabel('Time (days)')
ylabel('Energy loss (W/m^2)')
legend('Momentary energy loss')
hold off


ret = [(Time(1):Time(2):Time(3))', Energy, EnergyIn];

