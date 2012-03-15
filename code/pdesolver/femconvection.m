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

TneumannConditions = [10*h 0, NaN, 0;
		    -h, NaN, NaN, NaN];

TdirichletConditions = [NaN, NaN, 0, NaN];

UdirichletConditions = [NaN, NaN, NaN, 0];
UneumannConditions = [NaN, NaN, NaN, NaN;
		    NaN, NaN, NaN, NaN];

WdirichletConditions = [0, NaN, NaN, NaN];
WneumannConditions = [NaN, NaN, NaN, NaN;
                    NaN, NaN, NaN, NaN];



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
f = sparse(pCount,1);
GT = sparse(pCount, pCount);
Gu = sparse(pCount, pCount);
Gw = sparse(pCount, pCount);


divxx = sparse(pCount, pCount);
divxz = sparse(pCount, pCount);
divzx = sparse(pCount, pCount);
divzz = sparse(pCount, pCount);
area = sparse(1, tCount);

basedx = sparse(pCount, tCount);
basedz = sparse(pCount, tCount);

TensX = sptensor([pCount pCount pCount]);
TensZ = sptensor([pCount PCount pCount]);

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

neumann = e([1 2 5], find(~neunan));
neumannT = e([1 2 5], find(~neuTnan));
dirichlet = e([1 2 5], find(~dirnan));

%Create derivative matrix
fprintf(1, 'Calculating the derivatives of the basis...')

for k = 1:tCount
  basedx(t(1:3,k),k) = triderivative(p(corn), [1 2 3], [2 1 1], 1);
  basedz(t(1:3,k),k) = triderivative(p(corn), [1 2 3], [2 1 1], 2);
end

fprintf(1, ' done!\n');

%Calculate area of triangles

fprintf(1, 'Calculating the areas of the base triangles...');

for k = 1:tCount
  area = triarea(p(:,t(1:3,k)));
end

fprintf(1, ' done!\n');

%Assemble stiffness tensors
fprintf(1, 'Assembling stiffness tensors...');

for k = 1:tCount
  corn = t(1:3,k);
  
  vols = area(k)*(ones(3,3)/12 + eye(3,3)/12);
  %This formula has been taken from
  %http://users.wpi.edu/~sullivan/WebSite-ME515/Lectures/FiniteElement/Triangles/Triangles.htm
  
  dx = basedx(corn, k);
  dz = basedz(corn, k);
  
  TensX(corn(1),:,:) =  TensX(corn(1),:,:) +  vols*dx(1);
  TensX(corn(2),:,:) = TensX(corn(2),:,:) + vols*dx(2);
  TensX(corn(3),:,:) = TensX(corn(3),:,:) + vols*dx(3);

  TensZ(corn(1),:,:) = TensZ(corn(1),:,:) + vols*dz(1);
  TensZ(corn(2),:,:) = TensZ(corn(2),:,:) + vols*dz(2);
  TensZ(corn(3),:,:) = TensZ(corn(3),:,:) + vols*dz(3);
  
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
			    area(k)*basedx(t(1:3,k),t)*basedx(t(1:3,k), k)';

  divxz(t(1:3,k),t(1:3,k)) = divxz(t(1:3,k),t(1:3,k)) + ...
                            area(k)*basedx(t(1:3,k), k)*basedz(t(1:3,k),k)';

  divzz(t(1:3,k),t(1:3,k)) = divzz(t(1:3,k), t(1:3,k) + ...
			    area(k)*basedz(t(1:3,k))*basedz(t(1:3,k))';

end

divzx = divxz';

fprintf(1, ' done!\n');

%Assemble load vector

fprintf(1, 'Assembling load vector...');

for k = 1:tCount
  f(t(1:3,k),1) = f(t(1:3,k),1) + area(k)/3;
end

f = g*beta*f;

fprintf(1, ' done!\n');


%Enforce neumann conditions
fprintf(1, 'Enforcing neumann conditions...')


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

fprintf(1, ' done!\n')

%Move dirichlet nodes to RHS
fprintf(1, 'Enforcing dirichlet conditions...')
u(dirichlet(1,:)) = dirichletConditions(dirichlet(3,:));
u(dirichlet(2,:)) = dirichletConditions(dirichlet(3,:));

[int, ca, cb] = intersect(dirichlet(1,:), dirichlet(2,:));

ind = find(dirichlet(3,ca) ~= dirichlet(3,cb));

u(int(ind)) = 0.5*(dirichletConditions(dirichlet(3,ca(ind))) + ...
                   dirichletConditions(dirichlet(3,cb(ind))));

b = b - (A+Q) * u;
fprintf(1, ' done!\n');


%Solve system
fprintf(1, 'Solving system...');
FreeNodes=setdiff(1:pCount,unique(dirichlet));

u(FreeNodes) = (A(FreeNodes, FreeNodes)+Q(FreeNodes, FreeNodes))\ ...
    (b(FreeNodes) +G(FreeNodes));
fprintf(1, ' done!\n\n')
%Display
time = toc;

disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
disp(['Degrees of freedom: ' num2str(max(size(FreeNodes)))])
figure(1)
pdeplot(p,e,t,'xydata',u,'mesh',showgrid);

figure(2)
pdemesh(p,e,t)

figure(1)
