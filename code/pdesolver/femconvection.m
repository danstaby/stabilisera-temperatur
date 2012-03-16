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
border = [0,10,10,0;
    0,0,6,6];
h = 1e6;

g = 9.81;
beta = 3.67e-3; %### NOTE #### T-dependent. Fix later
alpha =1.9e-5; %Pressure and water dependant
rho = 1.2920; %Pressure and water dependant
nu = 13; %Temperature and pressure dependant
penalty = 1e6; %Which penalty?
Tref = 20;

TneumannConditions = [0, -1000000, 0, NaN];
TneumannTConditions = [NaN, NaN, NaN,NaN];

TdirichletConditions = [NaN, NaN, NaN, 20];

UdirichletConditions = [0, 0, 0, 0];
UneumannConditions = [NaN, NaN, NaN, NaN];
UneumannTConditions = [NaN, NaN, NaN, NaN];

WdirichletConditions = [0, 0, 0, 0];
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
LM = sparse(pCount, pCount);

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

Uneumann = e([1 2 5], find(~Uneunan));
UneumannT = e([1 2 5], find(~UneuTnan));
Udirichlet = e([1 2 5], find(~Udirnan));


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



voltemp = [2, 1, 1; 1,2,1;1,1,2]/12;

for k = 1:tCount
  %corn = t(1:3,k)';
  
  
  vols = full(area(k,1))*voltemp;
  %This formula has been taken from
  %http://users.wpi.edu/~sullivan/WebSite-ME515/Lectures/FiniteElement/Triangles/Triangles.htm
  
  
  for j = 1:3
    Tx(:,:,j) = sptensor(vols*basedx(t(j,k)', k));
    Tz(:,:,j) = sptensor(vols*basedz(t(j,k)', k));
  end
  

  TensX(t(1:3,k)', t(1:3,k)', t(1:3,k)') = TensX(t(1:3,k)', t(1:3,k)', t(1:3,k)') + Tx;
  TensZ(t(1:3,k)', t(1:3,k)',t(1:3,k)') = TensZ(t(1:3,k)', t(1:3,k)', t(1:3,k)') + Tz;
  
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
  
  fw(t(1:3,k),1) = fw(t(1:3,k),1) + -Tref*area(k)/3;
end

fw = g*beta*fw;

fprintf(1, ' done!\n');

%Assembling load matrix

fprintf(1, 'Assembling inhomogenous stiffness matrix...');

for j=1:tCount
  LM(t(1:3,j), t(1:3,j)) = LM(t(1:3,j), t(1:3,j)) - g*beta*area(j,1)*voltemp;
end

fprintf(1, ' done!\n');

%Enforce neumann conditions
fprintf(1, 'Enforcing neumann conditions...')


for j = 1 :  size(Tneumann,2)
  GT(Tneumann(1:2,j))=GT(Tneumann(1:2,j)) + ... 
      norm(p(:,Tneumann(1,j)) - ... 
	   p(:,Tneumann(2,j)))* ...
      ones(2,1)*TneumannConditions(Tneumann(3,j))/2;
end


for j = 1:size(Uneumann,2)
  Gu(Uneumann(1:2,j))=Gu(Uneumann(1:2,j)) + ... 
      norm(p(:,Uneumann(1,j)) - ... 
           p(:,Uneumann(2,j)))* ...
      ones(2,1)*UneumannConditions(Uneumann(3,j))/2;
end

for j = 1:size(Wneumann,2)
  Gw(Wneumann(1:2,j))=Gw(Wneumann(1:2,j)) + ... 
      norm(p(:,Wneumann(1,j)) - ... 
           p(:,Wneumann(2,j)))* ...
      ones(2,1)*WneumannConditions(Wneumann(3,j))/2;
end

for k = 1:size(TneumannT,2)
  L = norm(p(:,TneumannT(1,k)) - p(:,TneumannT(2,k)));
  
  QT(TneumannT([1 2],k), TneumannT([1 2],k)) = QT(TneumannT([1 2],k),TneumannT([1 2],k)) ...
      - TneumannTConditions(TneumannT(3,k))*L*[2,1;1,2]/6;
end

%QT

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
     squeeze(sptensor(nu*A*uu + penalty*(divxx*uu + divxz*wu)/rho));
fw = fw - ttv(ttv(TensX, wu,3),uu,2) - ttv(ttv(TensZ, wu, 3),wu,2) - ...
     (squeeze(sptensor(nu*A*wu + LM*Tu + penalty*(divzx*uu + divzz*wu)/rho)));

%b = b - (A+Q) * u;
fprintf(1, ' done!\n');


%Solve system
fprintf(1, 'Preparing for solution...');

ITERMAX = 40;
epsilon = 6e-5;

TFreeNodes = setdiff(1:pCount,unique(Tdirichlet));
uFreeNodes = setdiff(1:pCount, unique(Udirichlet));
wFreeNodes = setdiff(1:pCount, unique(Wdirichlet));

TFreeCount = size(TFreeNodes,2);
uFreeCount = size(uFreeNodes,2);
wFreeCount = size(wFreeNodes,2);

uStart = 1;
uEnd = uFreeCount;
wStart = uEnd + 1;
wEnd = wStart + wFreeCount - 1;
TStart = wEnd + 1;
TEnd = TStart + TFreeCount - 1;

resvec = zeros(uFreeCount+wFreeCount + TFreeCount, 1);
x = zeros(uFreeCount+wFreeCount+TFreeCount,1);
fprintf(1, ' done!\n');

fprintf(1, '\nSolving system with Newton-Raphson:\n');
disp(['Max iterations: ' num2str(ITERMAX) ', Epsilon: ' num2str(epsilon)])
fprintf(1, '\n')
res = epsilon + 1;
iterCount = 0;


order = 1;
lastres = res + 1;
%Make initial guess
uGuess = zeros(uFreeCount,1);
TGuess = Tref*ones(TFreeCount,1);
wGuess = zeros(wFreeCount,1);

fprintf(1,'ITERATION\tTRes\t\tuRes\t\twRes\n');

while((res > epsilon || iterCount < 3) && iterCount < ITERMAX)
  
  %Update function value
  
  resvec(uStart:uEnd,1) = double(ttv(ttv(TensX(uFreeNodes, uFreeNodes, uFreeNodes), uGuess,3),uGuess,2) + ...
      ttv(ttv(TensZ(uFreeNodes, wFreeNodes, uFreeNodes), uGuess, 3),wGuess,2)) +  ...
      nu*A(uFreeNodes, uFreeNodes)*uGuess + penalty*(divxx(uFreeNodes, uFreeNodes)*uGuess + ...
      divxz(uFreeNodes, wFreeNodes)*wGuess)/rho - ...
      double(fu(uFreeNodes));

  
  resvec(wStart:wEnd,1) = double(ttv(ttv(TensX(wFreeNodes, uFreeNodes,wFreeNodes), wGuess,3),uGuess,2) + ...
      ttv(ttv(TensZ(wFreeNodes, wFreeNodes, wFreeNodes), wGuess,3),wGuess,2)) + ... 
      nu*A(wFreeNodes, wFreeNodes)*wGuess + ...
      penalty*(divzx(wFreeNodes,uFreeNodes)*uGuess + ...
	       divzz(wFreeNodes,wFreeNodes)*wGuess)/rho + LM(wFreeNodes,TFreeNodes)*TGuess - double(fw(wFreeNodes));



						
  resvec(TStart:TEnd,1) = double(ttv(ttv(TensX(TFreeNodes, uFreeNodes,TFreeNodes),TGuess,3),uGuess,2)+...
      ttv(ttv(TensZ(TFreeNodes, wFreeNodes, TFreeNodes), TGuess,3), wGuess,2)) + ...
      alpha*((A(TFreeNodes, TFreeNodes)+QT(TFreeNodes,TFreeNodes))*TGuess) ...
	     -double(fT(TFreeNodes));

  %Calculate Jacobian

  Jacobian = sparse(uFreeCount + wFreeCount + TFreeCount,...
		    uFreeCount+wFreeCount + TFreeCount);
  
    
  %Jacobian(uStart:uEnd, uStart:uEnd) = double(ttv(TensX, uGuess),2)
  Jacobian(uStart:uEnd,uStart:uEnd) = double(ttv(TensZ(uFreeNodes, wFreeNodes,uFreeNodes),wGuess,2) ...
				    + ttv(TensX(uFreeNodes,uFreeNodes, uFreeNodes), uGuess, 2) ...
                                    + ttv(TensX(uFreeNodes,uFreeNodes,uFreeNodes), uGuess,3)) ...
                                    + penalty*divxx(uFreeNodes,uFreeNodes)/rho...
                                    + nu*A(uFreeNodes,uFreeNodes);

  Jacobian(uStart:uEnd,wStart:wEnd) = double(ttv(TensZ(uFreeNodes,wFreeNodes,uFreeNodes),uGuess,3))...
                                    + penalty*divxz(uFreeNodes,wFreeNodes)/rho;

  Jacobian(wStart:wEnd,wStart:wEnd) = double(ttv(TensX(wFreeNodes, uFreeNodes, wFreeNodes),uGuess,2) ...
				    + ttv(TensZ(wFreeNodes,wFreeNodes,wFreeNodes),wGuess,2) ...
		                    + ttv(TensZ(wFreeNodes, wFreeNodes, wFreeNodes), wGuess,3))...
                                    + penalty*divzz(wFreeNodes,wFreeNodes)/rho ...
                                    + nu*A(wFreeNodes, wFreeNodes);

  Jacobian(wStart:wEnd,uStart:uEnd) = double(ttv(TensX(wFreeNodes,uFreeNodes,wFreeNodes),wGuess,3))...
                                    + penalty*divzx(wFreeNodes,uFreeNodes)/rho;
  
  Jacobian(wStart:wEnd,TStart:TEnd) = LM(wFreeNodes,TFreeNodes);

  Jacobian(TStart:TEnd,TStart:TEnd) = double(ttv(TensX(TFreeNodes,uFreeNodes,TFreeNodes),uGuess,2)...
				    + ttv(TensZ(TFreeNodes,wFreeNodes,TFreeNodes),wGuess,2)) ...
                                    + alpha*(A(TFreeNodes,TFreeNodes)+QT(TFreeNodes, TFreeNodes));

  Jacobian(TStart:TEnd,uStart:uEnd) = double(ttv(TensX(TFreeNodes,uFreeNodes, uFreeNodes),uGuess,3));

  Jacobian(TStart:TEnd,wStart:wEnd) = double(ttv(TensZ(TFreeNodes,wFreeNodes,wFreeNodes),wGuess,3));

  %Note the factor 2 in the above expressions. That factor deals
  %with the square terms in the LHS of the equation

  %Solve for dx
    
  %[SU SS SV] = svds(uJacob);
  %dx(uStart:uEnd,1) = -SU*pinv(SS)*SV'*resvec(uStart:uEnd);
  %[SU SS SV] = svds(wJacob);
  %dx(wStart:wEnd,1) = -SU*pinv(SS)*SV'*resvec(wStart:wEnd);

  %[SU SS SV] = svds(TJacob);
  %dx(TStart:TEnd,1) = -SU*pinv(SS)*SV'*resvec(TStart:TEnd);

  %dx =  U*pinv(S)*V'*(-resvec); %Jacobian\(-resvec);
  %JacInv = inv(Jacobian);
  dx = Jacobian\(-resvec);

  %dx(uStart:uEnd,1) = -uJacob\resvec(uStart:uEnd);
  %dx(wStart:wEnd,1) = -wJacob\resvec(wStart:wEnd);
  %dx(TStart:TEnd,1) = -TJacob\resvec(TStart:TEnd);
  
  res = norm(resvec);

  %norm(wJacob)*norm(inv(wJacob))
  
  
  lastres = res;
  %Update vectors
  x = x + order*dx;
  uGuess = x(uStart:uEnd,1);
  wGuess = x(wStart:wEnd,1);
  TGuess = x(TStart:TEnd,1);
  res=norm(resvec);
  iterCount = iterCount + 1;
  

  fprintf(1, '%i\t\t%6.4g\t%6.4g\t%6.4g\n', iterCount, norm(resvec(TStart:TEnd)),...
	  norm(resvec(uStart:uEnd)), norm(resvec(wStart:wEnd)))
  %disp(['Iteration ' num2str(iterCount) ' done with residual ' ...
	%num2str(res)])
end

if(res < epsilon)
  fprintf(1, '\nFound convergence\n\n');
else
  fprintf(1, '\nNewton Raphson failed!\n\n');
end

Tu(TFreeNodes) = TGuess;
uu(uFreeNodes) = uGuess;
wu(wFreeNodes) = wGuess;
%u(FreeNodes) = (A(FreeNodes, FreeNodes)+Q(FreeNodes, FreeNodes))\ ...
%    (b(FreeNodes) +G(FreeNodes));

%Display
time = toc;

disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
disp(['Degrees of freedom: ' num2str(TFreeCount + wFreeCount + uFreeCount)])

figure(3)
pdemesh(p,e,t)

figure(2)

quiver(p(1,:)', p(2,:)', uu, wu);
figure(1)
%pdeplot(p,e,t,'xydata',Tu,'mesh',showgrid);
tricontour(p',t(1:3,:)', Tu, 10)

figure(4)
pdeplot(p,e,t,'xydata',Tu,'mesh',showgrid);
