function ret = femconvection(refinements, showgrid)
%Proof of concept application for heat equation in one dimension

if(nargin < 2)
  showgrid = 'off'
end

if(nargin < 1)
  refinements = 0;
end

tic;

dims = 1;
border = [0,1,1,0;
    0,0,1,1];
h = 1;

g = 9.81;
beta = 3.67e-3; %### NOTE #### T-dependent. Fix later
alpha = 1.9e-5; %Pressure and water dependant
rho = 1.2920; %Pressure and water dependant
nu = 13; %Temperature and pressure dependant
penalty = 1e7; %Which penalty?


TneumannConditions = [10*h 0, NaN, 0];
TneumannTConditions = [-h, NaN, NaN, NaN];

TdirichletConditions = [NaN, NaN, 0, NaN];

UdirichletConditions = [NaN, NaN, NaN, 0];
UneumannConditions = [NaN, NaN, NaN, NaN];
UneumannTConditions = [NaN, NaN, NaN, NaN];

WdirichletConditions = [0, NaN, NaN, NaN];
WneumannConditions = [NaN, NaN, NaN, NaN];

WneumannTConditions = [NaN, NaN, NaN, NaN];


%dbound = {{'20'}, {NaN}, {NaN}, {'0'}};
%nbound = {{NaN}, {NaN}, {'0'}, {NaN}};

fprintf(1,'Creating mesh...')
gd = [2,max(size(border)), border(1,:), border(2,:)]';
dl = decsg(gd);
[p, e, t] = initmesh(dl);
for n = 1:refinements
  [p, e, t] = refinemesh(dl, p, e, t);
end

%bMatrix = createBoundaryMatrix(dims,max(size(border)), dbound, nbound)
%u = assempde(bMatrix, p, e, t, 1, 0, 0);
%figure(3)
%pdeplot(p,e,t,'xydata',u,'mesh','off');

%Initiate variables
tCount = size(t,2);
pCount = size(p,2);
eCount = size(e,3);
A = sparse(pCount, pCount);
fu = sparse(pCount,1);
fw = sparse(pCount, 1);
uu = sparse(pCount, 1);
wu = sparse(pCount, 1);
Tu = sparse(pCount, 1);
fT = sparse(pCount, 1);
GT = sparse(pCount, 1);
Gu = sparse(pCount, 1);
Gw = sparse(pCount, 1);
QT = sparse(pCount, pCount);

divxx = sparse(pCount, pCount);
divxz = sparse(pCount, pCount);
divzx = sparse(pCount, pCount);
divzz = sparse(pCount, pCount);
area = sparse(1, tCount);

basedx = sparse(pCount, tCount);
basedz = sparse(pCount, tCount);

TensX = sptensor([pCount pCount pCount]);
TensZ = sptensor([pCount pCount pCount]);

fprintf(1,' done!\n');

Tneunan = isnan(TneumannConditions(e(5,:)));
TneuTnan = isnan(TneumannTConditions(e(5,:)));
Tdirnan = isnan(TdirichletConditions(e(5,:)));

Uneunan = isnan(UneumannConditions(e(5,:)));
UneuTnan = isnan(UneumannTConditions(e(5,:)));
Udirnan = isnan(UdirichletConditions(e(5,:)));

Wneunan = isnan(WneumannConditions(e(5,:)));
WneuTnan = isnan(WneumannTConditions(e(5,:)));
Wdirnan = isnan(WdirichletConditions(e(5,:)));

Tneumann = e([1 2 5], find(~Tneunan));
TneumannT = e([1 2 5], find(~TneuTnan));
Tdirichlet = e([1 2 5], find(~Tdirnan));

Wneumann = e([1 2 5], find(~Wneunan));
WneumannT = e([1 2 5], find(~WneuTnan));
Wdirichlet = e([1 2 5], find(~Wdirnan));

Uneumann = e([1 2 5], find(~Wneunan));
UneumannT = e([1 2 5], find(~WneuTnan));
Udirichlet = e([1 2 5], find(~Wdirnan));


%Create derivative matrix
fprintf(1, 'Calculating the derivatives of the basis...')

for k = 1:tCount
  deriv = triderivative(p(:,t(1:3,k)));
  basedx(t(1:3,k),k) = deriv(:,1);
  basedz(t(1:3,k),k) = deriv(:,2);
end

fprintf(1, ' done!\n');

%Calculate area of triangles

fprintf(1, 'Calculating the areas of the base triangles...');

for k = 1:tCount
  area(k,1) = triarea(p(:,t(1:3,k)));
end

fprintf(1, ' done!\n');

%Assemble stiffness tensors
fprintf(1, 'Assembling stiffness tensors...');

Tx = sptensor(3,3,3);
Tz = sptensor(3,3,3);

for k = 1:tCount
  corn = t(1:3,k)';
  
  vols = sptensor(area(k,1)*(ones(3,3)/12 + eye(3,3)/12));
  %This formula has been taken from
  %http://users.wpi.edu/~sullivan/WebSite-ME515/Lectures/FiniteElement/Triangles/Triangles.htm
  
  
  for j = 1:3
    Tx(:,:,j) = vols*basedx(corn(j), k);
    Tz(:,:,j) = vols*basedz(corn(j), k);
  end

  TensX(corn, corn, corn) = TensX(corn, corn, corn) + Tx;
  TensZ(corn, corn, corn) = TensZ(corn, corn, corn) + Tz;
  
end
fprintf(1, ' done!\n');

%Assemble stiffness matrix
fprintf(1,'Assembling stiffness matrix...');
for k = 1:tCount
  A(t(1:3,k), t(1:3,k)) = A(t(1:3,k),t(1:3,k)) + laplacestiff(p(:, t(1:3,k)));
end
fprintf(1, ' done!\n')

%Assemble divergence stiffness matrices

fprintf(1, 'Assembling divergence stiffness matrices...');

for k = 1:tCount
  divxx(t(1:3,k),t(1:3,k)) = divxx(t(1:3,k),t(1:3,k)) + ...
			    area(k)*basedx(t(1:3,k), k)*basedx(t(1:3,k), k)';

  divxz(t(1:3,k),t(1:3,k)) = divxz(t(1:3,k),t(1:3,k)) + ...
                            area(k)*basedx(t(1:3,k), k)*basedz(t(1:3,k),k)';

  divzz(t(1:3,k),t(1:3,k)) = divzz(t(1:3,k), t(1:3,k)) + ...
			    area(k)*basedz(t(1:3,k),k)*basedz(t(1:3,k),k)';

end

divzx = divxz';

fprintf(1, ' done!\n');

%Assemble load vector

fprintf(1, 'Assembling load vector...');

for k = 1:tCount
  
  fw(t(1:3,k),1) = fw(t(1:3,k),1) + area(k)/3;
end

fw = g*beta*fw;

fprintf(1, ' done!\n');


%Enforce neumann conditions
fprintf(1, 'Enforcing neumann conditions...')


for j = 1 :  size(Tneumann,2)
  GT(Tneumann(1:2,j))=GT(Tneumann(1:2,j)) + ... 
      norm(p(:,Tneumann(1,j)) - ... 
	   p(:,Tneumann(2,j))).* ...
      ones(2,1)*TneumannConditions(Tneumann(3,j))/2;
end


for j = 1:size(Uneumann,2)
  Gu(Uneumann(1:2,j))=Gu(Uneumann(1:2,j)) + ... 
      norm(p(:,Uneumann(1,j)) - ... 
           p(:,Uneumann(2,j))).* ...
      ones(2,1)*UneumannConditions(Uneumann(3,j))/2;
end

for j = 1:size(Wneumann,2)
  Gw(Wneumann(1:2,j))=Gw(Wneumann(1:2,j)) + ... 
      norm(p(:,Wneumann(1,j)) - ... 
           p(:,Wneumann(2,j))).* ...
      ones(2,1)*WneumannConditions(Wneumann(3,j))/2;
end

for k = 1:size(TneumannT,2)
  L = norm(p(:,TneumannT(1,k)) - p(:,TneumannT(2,k)));
  
  QT(TneumannT([1 2],k), TneumannT([1 2],k)) = QT(TneumannT([1 2],k),TneumannT([1 2],k)) ...
      - TneumannTConditions(TneumannT(3,k))*L*[2,1;1,2]/6;
end

%Note that we don't have any speed dependance in our Neumann
%boundary conditions for the velocity vector.

%Also note that we are setting the divergence of the velocity
%vector to zero on the border. 

fprintf(1, ' done!\n')

%Move dirichlet nodes to RHS
fprintf(1, 'Enforcing dirichlet conditions...')
Tu(Tdirichlet(1,:)) = TdirichletConditions(Tdirichlet(3,:));
Tu(Tdirichlet(2,:)) = TdirichletConditions(Tdirichlet(3,:));
uu(Udirichlet(1,:)) = UdirichletConditions(Udirichlet(3,:));
uu(Udirichlet(2,:)) = UdirichletConditions(Udirichlet(3,:));
wu(Wdirichlet(1,:)) = WdirichletConditions(Wdirichlet(3,:));
wu(Wdirichlet(2,:)) = WdirichletConditions(Wdirichlet(3,:));

[int, ca, cb] = intersect(Tdirichlet(1,:), Tdirichlet(2,:));

ind = find(Tdirichlet(3,ca) ~= Tdirichlet(3,cb));

Tu(int(ind)) = 0.5*(TdirichletConditions(Tdirichlet(3,ca(ind))) + ...
                   TdirichletConditions(Tdirichlet(3,cb(ind))));


[int, ca, cb] = intersect(Udirichlet(1,:), Udirichlet(2,:));

ind = find(Udirichlet(3,ca) ~= Udirichlet(3,cb));

uu(int(ind)) = 0.5*(UdirichletConditions(Udirichlet(3,ca(ind))) + ...
                   UdirichletConditions(Udirichlet(3,cb(ind))));

[int, ca, cb] = intersect(Wdirichlet(1,:), Wdirichlet(2,:));

ind = find(Wdirichlet(3,ca) ~= Wdirichlet(3,cb));

Wu(int(ind)) = 0.5*(WdirichletConditions(Wdirichlet(3,ca(ind))) + ...
                   WdirichletConditions(Wdirichlet(3,cb(ind))));



%Update load vectors with dirichlet conditions

fT = squeeze(sptensor(fT));
fu = squeeze(sptensor(fu));
fw = squeeze(sptensor(fw));


fT = fT - squeeze(sptensor(alpha*(A+QT)*Tu)) - ttv(ttv(TensX, Tu,3),uu,2) - ttv(ttv(TensZ, Tu, 3),wu,2);
fu = fu - ttv(ttv(TensX, uu,3),uu,2) - ttv(ttv(TensZ, uu, 3),wu,2) - ...
     squeeze(sptensor(nu*A*uu - divxx*uu - divxz*wu));
fw = fw - ttv(ttv(TensX, wu,3),uu,2) - ttv(ttv(TensZ, wu, 3),wu,2) - ...
     squeeze(sptensor(nu*A*wu - divzx*uu- divzz*wu));

%b = b - (A+Q) * u;
fprintf(1, ' done!\n');


%Solve system
fprintf(1, 'Solving system...');
%FreeNodes=setdiff(1:pCount,unique(dirichlet));

%u(FreeNodes) = (A(FreeNodes, FreeNodes)+Q(FreeNodes, FreeNodes))\ ...
%    (b(FreeNodes) +G(FreeNodes));
fprintf(1, ' done!\n\n')
%Display
time = toc;

disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
%disp(['Degrees of freedom: ' num2str(max(size(FreeNodes)))])
%figure(1)
%pdeplot(p,e,t,'xydata',u,'mesh',showgrid);

%figure(2)
pdemesh(p,e,t)

figure(1)
