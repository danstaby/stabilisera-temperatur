function [outData] = groundheatfemtransientanalys(refinements)
%outData = groundheatfemtransientanalys(refinements)
%
%Application to calculate the heat flow through the ground.
%refinements indicate how many refinemesh that shoud be run.

global tCoef

tCoef = getMeanTemp(0); %Make sine interpolation of average
                        %temperature data. tCoef = [a,b,c] with
			%T = a + b*cos(k*t) + c*sin(k*t) which is
                        %periodic for t = 365. 

if(nargin < 1)
  refinements = 0;
end

tic; %Start timer


kelvin = 273.15;

dims = 1;
border = [0,-40,-40, 52,52,12,  12,    0; %Create boundary.
	  0,  0,-30,-30, 0, 0,-0.5, -0.5]; %[X coordinates; Y coordinates]



h = 15.5; %Convective constant of the air
kGranite = 2.2;
alpha = kGranite/(790*2691); %Thermal diffusivity of granite
Tref = 10+kelvin;
Tin = 20+kelvin;


Tlow = -30+kelvin;
Thigh = 30 + kelvin;
Tstep = 0.5;

fprintf(1,'Creating mesh...')
gd = [2,max(size(border)), border(1,:), border(2,:)]';
dl = decsg(gd);

%Create mesh
[p, e, t] = initmesh(dl);
for n = 1:refinements
  [p, e, t] = refinemesh(dl, p, e, t);
end


%Initiate variables
tCount = size(t,2);
pCount = size(p,2);
eCount = size(e,2);
A = sparse(pCount, pCount);
M = sparse(pCount, pCount);
Q = sparse(pCount, pCount);
G = sparse(pCount, 1);
b = sparse(pCount, 1);
u = sparse(pCount, 1);
ab = zeros(pCount,1);
area = zeros(tCount,1);
lambda = zeros(pCount,1);
Const = zeros(pCount,1);
fprintf(1,' done!\n');

%Assemble stiffness matrix
fprintf(1,'Assembling stiffness matrix...');
for k = 1:tCount
  A(t(1:3,k), t(1:3,k)) = A(t(1:3,k),t(1:3,k)) ...
                        + alpha*laplacestiff(p(:, t(1:3,k)));
end

%Calculate the area of the triangles
for k = 1:tCount
  area(k,1) = triarea(p(:,t(1:3,k)));
end

%Prepare the mass matrix
for k = 1:tCount
  M(t(1:3,k), t(1:3,k)) = M(t(1:3,k), t(1:3,k)) ...
                        + area(k)*[2,1,1;1,2,1;1,1,2]/12;
end

fprintf(1, ' done!\n')

n = 0;
Tmean = zeros(3,floor((Thigh-Tlow)/Tstep));
pts = p(1:2, e(1,:)) - p(1:2, e(2,:));
LengthBound = sqrt(diag(pts'*pts));

Tsave = 0 + kelvin;


Tref = 20 + kelvin; %The temperature inside the building

Uconc = 0.7/(0.45); %The U-value of the foundation


uLast = (9+kelvin)*sparse(ones(pCount,1)); %Set initial values


%Settings vetor for the dirichlet conditions

dirichletConditions = [NaN, NaN, NaN,    NaN,    NaN,    NaN, NaN, NaN]; 
dirnan = isnan(dirichletConditions(e(5,:)));
dirichlet = e([1 2 5], find(~dirnan));
Free=setdiff(1:pCount,unique(dirichlet));


fprintf(1, 'Inversing matrices... ')
Minv(Free,Free) = inv(M(Free,Free));
MinvA(Free,Free) = Minv(Free,Free)*A(Free,Free);
fprintf(1, 'done!\n')


Q = sparse(pCount, pCount);
g = sparse(pCount, 1);
b = sparse(pCount, 1);
u = sparse(pCount, 1);





%Settings vectors for the neumann conditions.
% The coefficients are
% dT/dn = neumannConditions + T*neumannTConditions +
%neumannCosConditions*cos(omega*t) + neumannSinConditions*sin(omega*t)

neumannConditions =   [0  ,   0,   0, Uconc*Tref, Uconc*Tref, Uconc*Tref, ...
		       h*(tCoef(1)+kelvin), h*(tCoef(1)+kelvin)]/kGranite;
neumannTConditions =  [0, 0, 0,-Uconc,-Uconc,-Uconc,-h,-h]/kGranite;

neumannSinConditions = [NaN,NaN,NaN,NaN,NaN,NaN, ...
		        h*tCoef(3),h*tCoef(3)]/kGranite;
neumannCosConditions = [NaN,NaN,NaN,NaN,NaN,NaN,...
		        h*tCoef(2),h*tCoef(2)]/kGranite;


period = [0,0,0,0,0,0,365*24*3600, 365*24*3600];
omega = 2*pi/(365*24*3600);


neunan = isnan(neumannConditions(e(5,:)));
neuTnan = isnan(neumannTConditions(e(5,:)));


neumann = e([1 2 5], find(~neunan));
neumannT = e([1 2 5], find(~neuTnan));



%Enforce neumann conditions

for j = 1 :  size(neumann,2)
  g(neumann(1:2,j))=g(neumann(1:2,j)) + ... 
      norm(p(:,neumann(1,j)) - ... 
	   p(:,neumann(2,j))).* ...
      alpha*ones(2,1)*neumannConditions(neumann(3,j))/2;
end


for k = 1:size(neumannT,2)
  L = norm(p(:,neumannT(1,k)) - p(:,neumannT(2,k)));
  
  Q(neumannT([1 2],k), neumannT([1 2],k)) = ...
      Q(neumannT([1 2],k),neumannT([1 2],k)) ...
      - alpha*neumannTConditions(neumannT(3,k))*L*[2,1;1,2]/6;
end



%Move dirichlet nodes to RHS

u(dirichlet(1,:)) = dirichletConditions(dirichlet(3,:));
u(dirichlet(2,:)) = dirichletConditions(dirichlet(3,:));

[int, ca, cb] = intersect(dirichlet(1,:), dirichlet(2,:));

ind = find(dirichlet(3,ca) ~= dirichlet(3,cb));

u(int(ind)) = 0.5*(dirichletConditions(dirichlet(3,ca(ind))) + ...
		   dirichletConditions(dirichlet(3,cb(ind))));

b = b - (A+Q) * u + g;



%Solve system
fprintf(1,'Creating eigen value decomoposition... ')

[VecTemp lambdaTemp] = eig(full(MinvA(Free,Free)+Minv(Free,Free)* ...
				Q(Free,Free))); %Calculate
                                                %eigenvectors and
                                                %eigennumbers of
                                                %the problem


Vec(Free,Free) = VecTemp;
lambda(Free) = diag(lambdaTemp);
Vinv = zeros(pCount,pCount);
Vinv(Free,Free) = inv(Vec(Free,Free));
fprintf(1, ' done!\n')

eConst = sparse(pCount,1);
fConst = sparse(pCount,1);


neunanS = isnan(neumannSinConditions(e(5,:)));
neunanC = isnan(neumannCosConditions(e(5,:)));
neumannS = e([1 2 5], find(~neunanS));
neumannC = e([1 2 5], find(~neunanC));


%Enforce the sin and cos terms in the neumann conditions
for j = 1 :  size(neumannC,2)
  eConst(neumannC(1:2,j))=eConst(neumannC(1:2,j)) + ... 
      norm(p(:,neumannC(1,j)) - ... 
           p(:,neumannC(2,j))).* ...
      alpha*ones(2,1)*neumannCosConditions(neumannC(3,j))/2;
end

for j = 1 :  size(neumannC,2)
  fConst(neumannS(1:2,j))=fConst(neumannS(1:2,j)) + ... 
      norm(p(:,neumannS(1,j)) - ... 
           p(:,neumannS(2,j))).* ...
      alpha*ones(2,1)*neumannSinConditions(neumannS(3,j))/2;
end



%We will calculate a method of lines solution on the form
%T = Vec*(Y_h + Y_p)
%with
%Y_h =  Const*exp(-lambda*T)
%Y_p = aConst+bConst*cos(omega*t)+cConst*sin(omega*t)



%We have the MOL differential equation:
%T' + lambda*T = cAlpha + cBeta*cos(omega*t) + cGamma*sin(omega*t)

fprintf(1, 'Evaluating method of lines ODE...')

cAlpha = Vinv*Minv*b; 
cBeta = Vinv*Minv*eConst;
cGamma = Vinv*Minv*fConst;

%Calculate the coefficients in the particular solution.
aConst = cAlpha./lambda;
bConst = (lambda.*cBeta-omega*cGamma)./(omega^2+lambda.^2);
cConst = (cGamma+omega*bConst)./lambda;


Const = Vinv*uLast-aConst-bConst; %Set the constant before the
                                  %exponential term to match the
                                  %initial conditions.

fprintf(1, ' done!\n')

Tmean = zeros(3,1);
outData = [0:(3600*24):365*3600*24]'; %Initiate return data vectors
outData = [outData, zeros(size(outData,1),1)];

Lengths = zeros(3,1);
boundoffs = 3; %This is boundary dependant. Change if you're using
               %a new boundary.
for eid = 1:3
    ind = find(~(e(5,:)-eid-boundoffs)); %Calculate the lengths of
    Lengts(eid) = sum(LengthBound(ind)); %the foundation boundaries
end

n = 0;

fprintf(1, 'Iterating time steps... ')
%Calculate the solution for all each day in a year
for tNow = 0:(3600*24):365*3600*24
  n = n +1;
  u(Free) = Vec*(aConst+bConst*cos(omega*tNow)+ cConst*sin(omega*tNow));
                        %Note that the exponential term has been
                        %set to zero which is the limit when t->Inf
  for eid = 1:3
    ind = find(~(e(5,:)-eid-boundoffs));
    Tmean(eid) = mean(u(ind)); %Calculate the mean temperature
  end

  
  outData(n,2) = -Uconc*(Lengts*Tmean/sum(Lengts)-Tin); %Calculate
							%the
							%energy loss
end

fprintf(1, 'done!\n')


%Display

time = toc; %Stop timer

disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
disp(['Degrees of freedom: ' num2str(max(size(Free)))])


figure(1)
plot(outData(:,1)/(24*3600), outData(:,2))
xlabel('Tid (dygn)')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 365])

figure(2)
tricontourf(p(1:2,:), t(1:3,:), Vec*(aConst+bConst)-kelvin, ...
	    [1 17]);
xlabel('Position(m)')
ylabel('Position(m)')
title('Temperatur (C)')
figure(3)
tricontourf(p(1:2,:), t(1:3,:), Vec*(aConst-bConst)-kelvin, ...
	    [1 17]);
xlabel('Position(m)')
ylabel('Position(m)')
title('Temperatur (C)')
