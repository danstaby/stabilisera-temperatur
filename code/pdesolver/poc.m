function ret = poc(refinements, showgrid)
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
neumannConditions = [10*h 0, NaN, 0];
neumannTConditions = [-h, NaN, NaN, NaN];
dirichletConditions = [NaN, NaN, 0, NaN];

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
Q = sparse(pCount, pCount);
G = sparse(pCount, 1);
b = sparse(pCount, 1);
u = sparse(pCount, 1);

fprintf(1,' done!\n');

neunan = isnan(neumannConditions(e(5,:)));
neuTnan = isnan(neumannTConditions(e(5,:)));
dirnan = isnan(dirichletConditions(e(5,:)));

neumann = e([1 2 5], find(~neunan));
neumannT = e([1 2 5], find(~neuTnan));
dirichlet = e([1 2 5], find(~dirnan));


%Assemble stiffness matrix
fprintf(1,'Assembling stiffness matrix...');
for k = 1:tCount
  A(t(1:3,k), t(1:3,k)) = A(t(1:3,k),t(1:3,k)) + laplacestiff(p(:, t(1:3,k)));
end
fprintf(1, ' done!\n')


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
