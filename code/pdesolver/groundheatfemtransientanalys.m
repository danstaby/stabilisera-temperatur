function [outData] = groundheatfemtransientanalys(refinements, showgrid, displayOff)
%groundheatfem(refinements, showgrid)
%
%Application to calculate the heat flow through the ground.
%refinements indicate how many refinemesh that shoud be run.

global tCoef

tCoef = getMeanTemp(0);

if(nargin < 2)
  showgrid = 'off'
end

if(nargin < 1)
  refinements = 0;
end

tic;

kelvin = 273.15;

dims = 1;
border = [0,-40,-40, 52,52,12,  12,    0;
	  0,  0,-30,-30, 0, 0,-0.5, -0.5];

PrepareTempInterpolation([6, 0;16,10]);

h = 15.5;
Ug = 2.4;
kGranite = 2.2;
alpha = kGranite/(790*2691);
Tref = 10+kelvin;
Tin = 20+kelvin;


Tlow = -30+kelvin;
Thigh = 30 + kelvin;
Tstep = 0.5;

fprintf(1,'Creating mesh...')
gd = [2,max(size(border)), border(1,:), border(2,:)]';
dl = decsg(gd);

%if(nargin == 2)
  [p, e, t] = initmesh(dl);
  for n = 1:refinements
    [p, e, t] = refinemesh(dl, p, e, t);
  end
%end
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
  A(t(1:3,k), t(1:3,k)) = A(t(1:3,k),t(1:3,k)) + alpha*laplacestiff(p(:, ...
						  t(1:3,k)));
end

for k = 1:tCount
  area(k,1) = triarea(p(:,t(1:3,k)));
end


for k = 1:tCount
  M(t(1:3,k), t(1:3,k)) = M(t(1:3,k), t(1:3,k)) + area(k)*[2,1,1;1,2,1;1,1,2]/12;
end

fprintf(1, ' done!\n')

n = 0;
Tmean = zeros(3,floor((Thigh-Tlow)/Tstep));

pts = p(1:2, e(1,:)) - p(1:2, e(2,:));
LengthBound = sqrt(diag(pts'*pts));

Tsave = 0 + kelvin;

Tref = 20 + kelvin;
Uconc = 0.7/(0.45);

%if(nargin == 2)
  uLast = (9+kelvin)*sparse(ones(pCount,1));
%else
%  uLast = uSave;
%end
dirichletConditions = [NaN, NaN, NaN,    NaN,    NaN,    NaN, NaN, NaN];
dirnan = isnan(dirichletConditions(e(5,:)));
dirichlet = e([1 2 5], find(~dirnan));
Free=setdiff(1:pCount,unique(dirichlet));

Minv(Free,Free) = inv(M(Free,Free));
MinvA(Free,Free) = inv(M(Free,Free))*A(Free,Free);

WallID = 2;

n = 0;
dt = 5*24*3600;

outData = [(5:5:10*365)'];
outData = [outData, zeros(size(outData,1),1)];
outData = [outData, zeros(size(outData,1),1)];

Q = sparse(pCount, pCount);
g = sparse(pCount, 1);
b = sparse(pCount, 1);
u = sparse(pCount, 1);

boundoffs = 3;
Lengts = zeros(3,1);
for eid = 1:3
  ind = find(~(e(5,:)-eid-boundoffs));

  Lengts(eid) = sum(LengthBound(ind));
end



  
u(:) = 0;
b(:) = 0;
Q(:,:) = 0;
g(:) = 0;

n = n + 1;
%To =Tout(tNow) + kelvin;

neumannConditions =   [0  ,   0,   0, Uconc*Tref, Uconc*Tref, Uconc*Tref, ...
		    h*(tCoef(1)+kelvin), h*(tCoef(1)+kelvin)]/kGranite;
neumannTConditions =  [0, 0, 0,     -Uconc,     -Uconc,     -Uconc,     -h,  -h]/kGranite;

neumannSinConditions = [NaN,NaN,NaN,NaN,NaN,NaN,h*tCoef(3),h*tCoef(3)]/kGranite;
neumannCosConditions = [NaN,NaN,NaN,NaN,NaN,NaN,h*tCoef(2),h*tCoef(2)]/kGranite;
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
  
  Q(neumannT([1 2],k), neumannT([1 2],k)) = Q(neumannT([1 2],k),neumannT([1 2],k)) ...
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

[VecTemp lambdaTemp] = eig(full(MinvA(Free,Free)+Minv(Free,Free)* ...
				Q(Free,Free)));
Vec(Free,Free) = VecTemp;
lambda(Free) = diag(lambdaTemp);
Vinv = zeros(pCount,pCount);
Vinv(Free,Free) = inv(Vec(Free,Free));




Ce = sparse(pCount,1);
Cf = sparse(pCount,1);
eConst = sparse(pCount,1);
fConst = sparse(pCount,1);

neunanS = isnan(neumannSinConditions(e(5,:)));
neunanC = isnan(neumannCosConditions(e(5,:)));


neumannS = e([1 2 5], find(~neunanS));
neumannC = e([1 2 5], find(~neunanC));


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





%[p(1:2,e(1,:)); e(5,:)]

cAlpha = Vinv*Minv*b;
cBeta = Vinv*Minv*eConst;
cGamma = Vinv*Minv*fConst;

aConst = cAlpha./lambda;

bConst = (lambda.*cBeta-omega*cGamma)./(omega^2+lambda.^2);
cConst = (cGamma+omega*bConst)./lambda;



%Const(Free) = inv(Vec(Free,Free))*((uLast(Free)-ab(Free)));

Const = Vinv*uLast-aConst-bConst;
tNow = 0.5*pi/omega;


u(Free) = Vec*(Const.*exp(-lambda*tNow) +...
	       aConst+bConst*cos(omega*tNow)+ cConst*sin(omega*tNow));

%Calculate mean temp and save data.
  
%outData(n,2) = mean(u(find(~(e(5,:)-WallID))));
%uLast = u;
%fprintf(1, '\r%5g%%    ', 100*tNow/(365*10));

%Tmean = zeros(3,1);
%for eid = 1:3
%  ind = find(~(e(5,:)-eid-boundoffs));
%  
%  Tmean(eid) = mean(u(ind));
%end

%outData(n,3) = -Uconc*(Lengts'*Tmean/sum(Lengts)-Tin);


%Lengts = zeros(3,1);



%for eid = 1:3
%    ind = find(~(e(5,:)-eid-boundoffs));
%    
%    Lengts(eid) = sum(LengthBound(ind));
%end

%TotalMean = Lengts'*Tmean/sum(Lengts);

%Display
time = toc;

disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
disp(['Degrees of freedom: ' num2str(max(size(Free)))])

if(displayOff == 0)
  %figure(1)
  %pdeplot(p,e,t,'xydata',u-kelvin,'mesh',showgrid);
  %hold off
  %plot(outData(:,1)/(365), outData(:,2)-kelvin)
  %xlabel('Tid (år)')
  %ylabel('Medeltemperatur (C)')
  

  %figure(2)
  %hold off
  %plot(outData(:,1)/(365), outData(:,3))
  %xlabel('Tid (År)')
  %ylabel('Kyleffekt (W m^{-2})')
  figure(3)
  hold off
  
  
  tricontourf(p(1:2,:),t(1:3,:), Vec*(aConst + bConst)-kelvin);
  
end

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


function ret = Tout(Time)
global tCoef

ret = tCoef(1)+tCoef(2)*cos(2*pi*Time/365)+tCoef(3)*sin(2*pi*Time/365);

function ret = ToutDaily(Time)
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
