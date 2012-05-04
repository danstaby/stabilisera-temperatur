function ret = groundheatfemiter(refinements, showgrid)
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
h = 15.5;
Ug = 2.4;
kGranite = 2.2;
Tref = 10+kelvin;
Tin = 20+kelvin;
Tout = 0+kelvin;

Tlow = -30+kelvin;
Thigh = 30 + kelvin;
Tstep = 0.5;


fprintf(1,'Creating mesh...')
gd = [2,max(size(border)), border(1,:), border(2,:)]';
dl = decsg(gd);
[p, e, t] = initmesh(dl);
for n = 1:refinements
  [p, e, t] = refinemesh(dl, p, e, t);
end

%Initiate variables
tCount = size(t,2);
pCount = size(p,2);
eCount = size(e,2);
A = sparse(pCount, pCount);
Q = sparse(pCount, pCount);
G = sparse(pCount, 1);
b = sparse(pCount, 1);
u = sparse(pCount, 1);

fprintf(1,' done!\n');

%Assemble stiffness matrix
fprintf(1,'Assembling stiffness matrix...');
for k = 1:tCount
  A(t(1:3,k), t(1:3,k)) = A(t(1:3,k),t(1:3,k)) + laplacestiff(p(:, ...
						  t(1:3,k)));
end
fprintf(1, ' done!\n')

n = 0;
Tmean = zeros(3,floor((Thigh-Tlow)/Tstep));

pts = p(1:2, e(1,:)) - p(1:2, e(2,:));
LengthBound = sqrt(diag(pts'*pts));

Tsave = 0 + kelvin;

for Tout = Tlow:Tstep:Thigh
  Q = sparse(pCount, pCount);
  G = sparse(pCount, 1);
  b = sparse(pCount, 1);
  u = sparse(pCount, 1);
  n = n + 1;
  Tref = mean([Tin, Tout]);
  neumannConditions =   [0  ,   NaN,   0, Ug*Tref, Ug*Tref, Ug*Tref, h*Tout, h*Tout]/kGranite;
  neumannTConditions =  [NaN, NaN, NaN,     -Ug,     -Ug,     -Ug,     -h,  -h]/kGranite;
  dirichletConditions = [NaN, kelvin + 8, NaN,    NaN,    NaN,    NaN, NaN, NaN];

  neunan = isnan(neumannConditions(e(5,:)));
  neuTnan = isnan(neumannTConditions(e(5,:)));
  dirnan = isnan(dirichletConditions(e(5,:)));

  neumann = e([1 2 5], find(~neunan));
  neumannT = e([1 2 5], find(~neuTnan));
  dirichlet = e([1 2 5], find(~dirnan));


  %Enforce neumann conditions
  

  for j = 1 :  size(neumann,2)
    G(neumann(1:2,j))=G(neumann(1:2,j)) + ... 
	norm(p(:,neumann(1,j)) - ... 
	     p(:,neumann(2,j))).* ...
	ones(2,1)*neumannConditions(neumann(3,j))/2;
  end


  for k = 1:size(neumannT,2)
    L = norm(p(:,neumannT(1,k)) - p(:,neumannT(2,k)));
  
    Q(neumannT([1 2],k), neumannT([1 2],k)) = Q(neumannT([1 2],k),neumannT([1 2],k)) ...
	- neumannTConditions(neumannT(3,k))*L*[2,1;1,2]/6;
  end

  

  %Move dirichlet nodes to RHS
  
  u(dirichlet(1,:)) = dirichletConditions(dirichlet(3,:));
  u(dirichlet(2,:)) = dirichletConditions(dirichlet(3,:));

  [int, ca, cb] = intersect(dirichlet(1,:), dirichlet(2,:));

  ind = find(dirichlet(3,ca) ~= dirichlet(3,cb));

  u(int(ind)) = 0.5*(dirichletConditions(dirichlet(3,ca(ind))) + ...
		     dirichletConditions(dirichlet(3,cb(ind))));

  b = b - (A+Q) * u;
  


  %Solve system
  
  FreeNodes=setdiff(1:pCount,unique(dirichlet));

  u(FreeNodes) = (A(FreeNodes, FreeNodes)+Q(FreeNodes, FreeNodes))\ ...
      (b(FreeNodes) +G(FreeNodes));
  
  
  boundoffs = 3;
  for eid = 1:3
    ind = find(~(e(5,:)-eid-boundoffs));
    
    
    Tmean(eid,n) = u(ind)'*LengthBound(ind)/sum(LengthBound(ind)); %Weighted mean
    
  end
  disp(['Tout = ' num2str(Tout) ', Tmean = ' num2str(Tmean(2,n))])
  if(Tout == Tsave)
    uSave = u;
  end
end

Lengts = zeros(3,1);



for eid = 1:3
    ind = find(~(e(5,:)-eid-boundoffs));
    
    Lengts(eid) = sum(LengthBound(ind));
end

TotalMean = Lengts'*Tmean/sum(Lengts);

%Display
time = toc;

disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
disp(['Degrees of freedom: ' num2str(max(size(FreeNodes)))])
figure(1)
%pdeplot(p,e,t,'xydata',u-kelvin,'mesh',showgrid);

plot((Tlow:Tstep:Thigh)-kelvin, Tmean(1,:)-kelvin)
hold on
plot((Tlow:Tstep:Thigh)-kelvin, Tmean(2,:)-kelvin, 'r')
plot((Tlow:Tstep:Thigh)-kelvin, Tmean(3,:)-kelvin, 'k')
hold off
xlabel('Temperature outside (C)')
ylabel('Mean temperature of boundary (C)')
legend('Right side', 'Lower side', 'Left side')

figure(2)
plot((Tlow:Tstep:Thigh)-kelvin, TotalMean-kelvin)
xlabel('Temperature ute (C)')
ylabel('Grundens medeltemperatur (C)')
xlim([-30, 30])
%figure(2)
%pdemesh(p,e,t)
%tricontourf(p(1:2,:),t(1:3,:),u-kelvin)

OurU = abs(Ug*(10+kelvin-TotalMean)./(TotalMean-(Tlow:Tstep:Thigh)));

figure(3)
plot((Tlow:Tstep:Thigh)-kelvin, OurU)
xlabel('Temperatur ute (C)')
ylabel('U av grunden')

OurU = OurU(1)

MeanLength = 2.2/OurU

%figure(1)
figure(4)
tricontourf(p(1:2,:),t(1:3,:),uSave-kelvin)
xlabel('Position (m)')
ylabel('Position (m)')
title('Temperatur (C)')