function h = tricontourf(p,t,u)
%h = tricountourf(p,t,u)
%
%This function plots a filled countour of the data in u aligned in
%the points p.
%
% p - row vectors of x and y coordinates for the points
% t - matrix of the triangulation
% u - the data that should be plotted


A = sparse(3*size(t,2),3*size(t,2));
b = zeros(3*size(t,2),1);
%Create linear spline interpolation

for k = 1:size(t,2)
  cor = t(:,k);
  A((3*k-2):(3*k),(3*k-2):(3*k)) = [1,1,1;p(:,cor)]';
  b((3*k-2):(3*k)) = u(cor);
end

spline = A\b;

%Create mesh

ymin = min(p(2,:));
xmin = min(p(1,:));
ymax = max(p(2,:));
xmax = max(p(1,:));

side = floor(sqrt(max(size(p))));

resolution = 2;

xvec = (xmax-xmin)*[0:resolution*side]/(resolution*side) + xmin;
yvec = (ymax-ymin)*[0:resolution*side]/(resolution*side) + ymin;

[X Y] = meshgrid(xvec,yvec);

%Put interpolation on the mesh
pCount = size(X,1)*size(X,2);
Z = zeros(1,pCount);

Z(:) = NaN; %Z data

P = [reshape(X,1,pCount); reshape(Y,1,pCount); zeros(1,pCount)];
%Point matrix
p = [p;zeros(1,size(p,2))];
fillMatrix = ones(1,pCount); %Matrix to convert vector to matrix

mask = zeros(3,pCount); %A mask to only get the z component of the
                        %cross product
mask(3,:) = 1;

for k=1:size(t,2)
  
  a = P-p(:,t(1,k))*fillMatrix; %Perform coordinate swap for easier compution
  b = (p(:,t(3,k))-p(:,t(1,k)))*fillMatrix;
  c = (p(:,t(2,k))-p(:,t(1,k)))*fillMatrix;
  
  cab = dot(cross(a,b),mask); %Perform triangulation
  cac = dot(cross(a,c),mask);
  cbc = dot(cross(b,c),mask);
  
  u = cac./cbc; %Create coordinate transformation
  v = -cab./cbc;
  
  id = find(u >= 0 & v >= 0 & (v+u) <= 1);
  
  if(isempty(id) == 0)
    Z(id) = dot([ones(1,size(id,2));P(1:2,id)],spline((3*k-2):(3*k))*fillMatrix(1,1:max(size(id))));
  end
end

Z = reshape(Z,size(X,1),size(X,2));

[Co Ho] = contourf(X,Y,Z);
%clabel(Co,Ho)
colorbar('location','southoutside')


