function ret = walltransientfem(Nodes)

Width = 0.5;
alpha = 5.2e-7;

Tend = 30000;
Tout = 10;
Tin = 20;
TStart = 0;
dt = 1;
L = Width/(Nodes-1);
xvec = [0:L:Width]';
uInitial = xvec*(Tin-TStart)/Width + TStart;
uLast = uInitial;

A = sparse(Nodes, Nodes);
M = sparse(Nodes, Nodes);
u = sparse(Nodes,1);
g = sparse(Nodes,1);
f = sparse(Nodes,1);
y = sparse(Nodes,1);
b = sparse(Nodes,1);
Vec = zeros(Nodes,Nodes);
lambda = zeros(Nodes,1);

dirbc = [1,Nodes;
	Tout,Tin];

%Assemble stiffness matrix

for n=1:(Nodes-1)
  A([n n+1], [n n+1]) = A([n n+1], [n n+1]) + [1/L,-1/L; -1/L, 1/L];
end
A = alpha*A;
%Assemble mass matrix

for n = 1:(Nodes-1)
  M([n n+1], [n n+1]) = M([n n+1], [n n+1]) + 4*[L/12, L/24; L/24, L/12];
end

Free = 2:(Nodes-1);
Const = sparse(Nodes,1);
Minv = sparse(Nodes,Nodes);
Vinv = sparse(Nodes,Nodes);
Ainv = sparse(Nodes,Nodes);

%mb = sparse(Nodes,1);
ab = sparse(Nodes,1);
figure(1)
plot(xvec, uInitial)
hold on

[VecTemp lambdaTemp] = eig(full(inv(M(Free,Free))*A(Free,Free)));
Vec(Free,Free) = VecTemp; 
lambda(Free) = diag(lambdaTemp);
Minv(Free,Free) = inv(M(Free,Free));
Vinv(Free,Free) = inv(Vec(Free,Free));

Ainv(Free,Free) = inv(A(Free,Free));
for t = dt:dt:Tend

  %Solve differential equations
  u(dirbc(1,:)) = dirbc(2,:)';
  b = f - A*u;
  %mb(:) = Minv*b;
  ab(Free) = Ainv(Free,Free)*b(Free);
  
  Const(Free) = Vinv(Free,Free)*((uLast(Free)-ab(Free))); %./lambda(Free))
  
  u(Free) = Vec(Free,Free)*(exp(-lambda(Free)*dt).*Const(Free)) + ab(Free);
  
  %u(Free) = uLast(Free) + dt*inv(M(Free,Free))*(b(Free)-A(Free,Free)*uLast(Free));
  
  uLast(:) = u;
  u(:) = 0;
  if(mod(t,1000) == 0)
    plot(xvec, uLast, 'r');
  end
  
  if(t == 1000)
    
    %plot(xvec,[10;Vec(Free,Free)*(exp(-lambda(Free)*30000).* ...
	%			 Const(Free)) + ab(Free);20], 'k')
	plot(xvec, uLast, 'm')
  end
end



plot(xvec,uLast, 'k');

hold off