function createMesh(fileName, refinements)
%createMesh(fileName, refinements)
%
%This function precalculates a new mesh and saves it to a file.
%It also calculates the stiffness matrices and the stiffness tensor
%to save time while solving a specified problem with Galerkin's method.
%
%fileName - The name of the file to save to
%refinements - How many times the mesh will be divided into smaller
%elements.


tic; %Start timer to time how long the execution of this script takes

border = [0,20,20,-6,0;
         0,0,22,22,18];

fprintf(1,'Creating mesh...')
gd = [2,max(size(border)), border(1,:), border(2,:)]';
dl = decsg(gd);
[p, e, t] = initmesh(dl); %Initialize mesh. Function from pde toolbox
for n = 1:refinements
  p = jigglemesh(p,e,t, 'Opt', 'off', 'Iter', 50); %Improve quality
                                                   %of the mesh

  [p, e, t] = refinemesh(dl, p, e, t); %Add more nodes and triangles
end

p = jigglemesh(p,e,t, 'Opt', 'off', 'Iter', 50);
fprintf(1, ' done\n');

tCount = size(t,2);
pCount = size(p,2);
eCount = size(e,2);

A = sparse(pCount, pCount);

divxx = sparse(pCount, pCount);
divxz = sparse(pCount, pCount);
divzx = sparse(pCount, pCount);
divzz = sparse(pCount, pCount);

area = sparse(tCount, 1);

basedx = sparse(pCount, tCount);
basedz = sparse(pCount, tCount);

TensX = sptensor([pCount pCount pCount]);
TensZ = sptensor([pCount pCount pCount]);

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
  TensX(t(1:3,k)', t(1:3,k)', t(1:3,k)') = TensX(t(1:3,k)', t(1:3,k)', ...
					       t(1:3,k)') + Tx;
  TensZ(t(1:3,k)', t(1:3,k)',t(1:3,k)') = TensZ(t(1:3,k)', t(1:3,k)', ...
						t(1:3,k)') + Tz;

end
fprintf(1, ' done!\n');

%Assemble stiffness matrix
fprintf(1,'Assembling stiffness matrix...');
for k = 1:tCount
  A(t(1:3,k), t(1:3,k)) = A(t(1:3,k),t(1:3,k)) + laplacestiff(p(:, ...
						  t(1:3,k)));
end
fprintf(1, ' done!\n')

%Assemble divergence stiffness matrices

fprintf(1, 'Assembling divergence stiffness matrices...');

for k = 1:tCount
  divxx(t(1:3,k),t(1:3,k)) = divxx(t(1:3,k),t(1:3,k)) + ...
                            area(k,1)*basedx(t(1:3,k), k)*basedx(t(1:3,k), ...
						  k)';

  divxz(t(1:3,k),t(1:3,k)) = divxz(t(1:3,k),t(1:3,k)) + ...
                            area(k,1)*basedx(t(1:3,k), k)*basedz(t(1: ...
						  3,k),k)';

  divzz(t(1:3,k),t(1:3,k)) = divzz(t(1:3,k), t(1:3,k)) + ...
                            area(k,1)*basedz(t(1:3,k),k)*basedz(t(1: ...
						  3,k),k)';

end

divzx = divxz';

fprintf(1, ' done!\n');

fprintf(1, 'Saving file... ');
save(fileName, 'A', 'area', 'basedx', 'basedz', 'TensX', 'TensZ', ...
     'divxx', 'divzx', 'divxz', 'divzz', 'eCount', 'pCount', 'tCount', ...
     'p', 'e', 't', 'voltemp');
fprintf(1, 'done!\n\n');

time = toc; %Stop timer


disp(['Execution time: ' num2str(time) ' s']);
disp(['Triangle count: ' num2str(tCount)])
disp(['Node count: ' num2str(pCount)])
