function ret = walltransientfem(Time, Nodes, Insulated)
% walltransientfem(Time, Nodes, Insulated)
%
% This Script calculates the transient temperatures in a wall with
% a specified amount of nodes. Insulated = 0 means uninsulated wall
% and Insulated = 1 mean that the wall is insulated.



if(nargin == 2)
  Insulated = 0
end


if( Insulated == 0)
  Material = [0.5, 5.2e-7, 0.6];


else
  Material = [0.5, 5.2e-7, 0.6;
	      0.1, 0.037/(18*840), 0.037];
  %18 kg/m^3
  %840 J/K kg
  %0.037 W/mK
end



Width = sum(Material(:,1));

%Put nodes on boundaries and on the borders between materials
NrNodes = floor(Nodes*Material(:,1)/Width)+1;
NrNodes(NrNodes < 2)  = 2;
MaterialCount = size(Material,1);
NodeCount = sum(NrNodes)+1;
Pts = zeros(NodeCount,1);
L = zeros(NodeCount,1);
uInitial = sparse(NodeCount,1);
alpha = zeros(NodeCount-1,1);
PrevEnd = 0;
PrevL = 0;

for n = 1:size(Material,1)
  Pts(PrevEnd+(1:NrNodes(n))) = PrevL + Material(n,1)*(0: ...
						  (NrNodes(n)-1))'/(NrNodes(n));
  alpha(PrevEnd+(1:NrNodes(n))) = Material(n,2);
  L(PrevEnd + (1:NrNodes(n))) = Material(n,1)/(NrNodes(n));
  PrevL = PrevL + Material(n,1);
  PrevEnd = PrevEnd + NrNodes(n);
end

Pts(sum(NrNodes)+1) = PrevL;

Tu = 10;
Tin = 20;
TStart = 0;



%Set initial value to the steady state solution
if(size(Material,1) > 1)
  [~, TempInit] = heatrod([TStart Tin], Material(:, [1 3]));
end
PrevEnd = 0;
Tlast = TStart;

for n = 1:(size(Material,1)-1)
  
  uInitial(PrevEnd+(1:NrNodes(n))) = (Pts(PrevEnd+(1:NrNodes(n)))-sum(L(1:(PrevEnd))))* ...
                                    (TempInit(n)-Tlast)/Material(n,1)+Tlast;
  Tlast = TempInit(n);
  PrevEnd = PrevEnd + NrNodes(n);
end

uInitial(PrevEnd+(1:(NrNodes(MaterialCount)))) = (Pts(PrevEnd+(1:(NrNodes(MaterialCount))))-...
					   sum(L(1:(PrevEnd))))*(Tin-Tlast)/Material(MaterialCount,1)...
	                                   +Tlast;
uInitial(NodeCount) = Tin;
uLast = uInitial;



A = sparse(NodeCount, NodeCount);
M = sparse(NodeCount, NodeCount);
u = sparse(NodeCount,1);
f = sparse(NodeCount,1);
b = sparse(Nodes,1);
Vec = zeros(Nodes,Nodes);
lambda = zeros(Nodes,1);


dirbc = [1,NodeCount;
	Tu,Tin];

%Assemble stiffness matrix

for n=1:(NodeCount-1)
  A([n n+1], [n n+1]) = A([n n+1], [n n+1]) + alpha(n)*[1/L(n),-1/L(n); -1/L(n), 1/L(n)];
end

%Assemble mass matrix

for n = 1:(NodeCount-1)
  M([n n+1], [n n+1]) = M([n n+1], [n n+1]) + L(n)*[1/3, 1/6; 1/6, 1/3];
end

Free = 2:(NodeCount-1);
Const = sparse(NodeCount,1);
Minv = sparse(NodeCount,NodeCount);
Vinv = sparse(NodeCount,NodeCount);
Ainv = sparse(NodeCount,NodeCount);

%mb = sparse(Nodes,1);
ab = sparse(NodeCount,1);
figure(1)
plot(Pts, uInitial)
hold on
[VecTemp lambdaTemp] = eig(full(inv(M(Free,Free))*A(Free,Free)));
Vec(Free,Free) = VecTemp; 
lambda(Free) = diag(lambdaTemp);
Minv(Free,Free) = inv(M(Free,Free));
Vinv(Free,Free) = inv(Vec(Free,Free));

Ainv(Free,Free) = inv(A(Free,Free));
for t = Time(1):Time(2):Time(3)

  %Solve differential equations
  u(dirbc(1,:)) = dirbc(2,:)';
  b = f - A*u;
  %mb(:) = Minv*b;
  ab(Free) = Ainv(Free,Free)*b(Free);
  
  Const(Free) = Vinv(Free,Free)*((uLast(Free)-ab(Free))); %./lambda(Free))
  
  u(Free) = Vec(Free,Free)*(exp(-lambda(Free)*Time(2)).*Const(Free)) + ab(Free);
  
  %u(Free) = uLast(Free) + dt*inv(M(Free,Free))*(b(Free)-A(Free,Free)*uLast(Free));
  
  uLast(:) = u;
  u(:) = 0;
  if(mod(t,1000) == 0)
    plot(Pts, uLast, 'r');
  end
  
  if(t == 1000)
    
    %plot(xvec,[10;Vec(Free,Free)*(exp(-lambda(Free)*30000).* ...
	%			 Const(Free)) + ab(Free);20], 'k')
	plot(Pts, uLast, 'm')
  end
end



plot(Pts,uLast, 'k');

hold off