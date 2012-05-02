function ret = femsolveriter(fileName, Tout, VarName, Range)
%femsolveriter(fileName, Tout, VarName, Range)
%
%This function solves the navier stokes equations for
%a specified problem area.
%fileName - path to file specifying problem area
%Tout - The temperature outside
%VarName - String that contains the name of the variable to varied
%Range - Range of the variable. [Start, Step, End]

global TensX TensZ divxx divzz divzx divxz A QT LM alpha beta nu rho penalty Gu Gw GT


h = 6.61;
g = 9.81;
beta = 3.67e-3; %### NOTE #### T-dependent. Fix later
alpha =1.9e-5; %Pressure and water dependant
rho = 1.2920; %Pressure and water dependant
nu = 1.3e-5;%1.3e-5;%1.3e-5; %Temperature and pressure dependant
penalty = 1e7; %Which penalty?
sigma = 5.670373e-8;
kelvin = 273.15;
Tref = Tout+kelvin;
T0 = Tref+3; %A guess for wall temperature
TL = -20+kelvin;
Tin = 20 + kelvin;
Uvalue = 1.1858; %U-Value of the south, east and west walls
kAir = 0.0243;

ITERMAX = 400; %Max iterations for Newton Raphsons method
epsilon = 1e-5; %Break point for NR-convergence 


WallID = 5; %Identifier of the edge that reresents the wall.

meantemp = 0;
temperr = 0.05;
newtemp = T0;

fprintf(1,'Loading problem ...');
 S = load(fileName); %Loading problem specific variables
 p = S.p;
 e = S.e;
 t = S.t;
 tCount = S.tCount;
 pCount = S.pCount;
 eCount = S.eCount;
 A = sparse(S.A);
 divxx = sparse(S.divxx);
 divxz = sparse(S.divxz);
 divzx = sparse(S.divzx);
 divzz = sparse(S.divzz);
 area = sparse(S.area);
 basedx = sparse(S.basedx);
 basedz = sparse(S.basedz);
 TensX = S.TensX;
 TensZ = S.TensZ;
 voltemp = S.voltemp;
 fprintf(1, ' done!\n');

 uu = sparse(pCount, 1); %Solution vectors
 wu = sparse(pCount, 1);
 Tu = sparse(pCount, 1);


AirH = zeros(size(Range,2),1);
iterN = 0;



for ourN = 1:size(Range,2)
  VarN = Range(ourN);
  
  iterN = iterN + 1;
  eval([VarName '=' num2str(VarN) ';']);
  %newtemp = meantemp;
  meantemp = T0 + temperr + 1;
  disp([VarName '=' num2str(VarN)])
  while(abs(meantemp - T0) > temperr) 
    T0 = newtemp;
    %bb = 50 + sigma*(0.8*(0.7*(Tref-30)^4 +  0.3*Tref^4) + 0.8*3*T0^4); %Q due to planck radiation
    %bbt = -0.8*4*sigma*T0^3; %T-dependance of Q due to planck radiation
			     %(first order tayolor polynominal)

			     %Boundary conditions. T = TdirichletConditions,
			     % dT/dn = TneumannConditions +
                             % T*TneumannTConditions
     
     Rair = sigma*Tref^4*(1-0.261*exp(-7.77e-4*(273-Tref)));
     Ramb = sigma*Tref^4;
     
     RwallT = -0.90*4*sigma*T0^3;
     Rwall = 0.9*3*sigma*T0^4; 

     TneumannConditions = [NaN, h*Tref, h*Tref, 0, Uvalue*Tin + 0.38*Rair+0.62*Ramb+Rwall]/kAir;
     TneumannTConditions = [NaN,-h, -h, NaN, -Uvalue + RwallT]/kAir;
     TdirichletConditions = [Tref, NaN, NaN, NaN, NaN];

     UdirichletConditions = [0, 0, 0, 0, 0];
     UneumannConditions = [NaN, NaN, NaN, NaN, NaN];
     UneumannTConditions = [NaN, NaN, NaN, NaN, NaN];

     WdirichletConditions = [0, 0, 0, 0 ,0 ];
     WneumannConditions = [NaN, NaN, NaN, NaN, NaN];
     WneumannTConditions = [NaN, NaN, NaN, NaN, NaN];

     fu = sparse(pCount,1); %Load vectors
     fw = sparse(pCount, 1);
     fT = sparse(pCount, 1);


     GT = sparse(pCount, 1); %Neumann vectors
     Gu = sparse(pCount, 1);
     Gw = sparse(pCount, 1);

     QT = sparse(pCount, pCount); %T dependant neumann matrix
     LM = sparse(pCount, pCount); %Inhomogenous stiffness matrix

     Tneunan = isnan(TneumannConditions(e(5,:))); %Find active boundary conditions
     TneuTnan = isnan(TneumannTConditions(e(5,:)));
     Tdirnan = isnan(TdirichletConditions(e(5,:)));

     Uneunan = isnan(UneumannConditions(e(5,:)));
     UneuTnan = isnan(UneumannTConditions(e(5,:)));
     Udirnan = isnan(UdirichletConditions(e(5,:)));

     Wneunan = isnan(WneumannConditions(e(5,:)));
     WneuTnan = isnan(WneumannTConditions(e(5,:)));
     Wdirnan = isnan(WdirichletConditions(e(5,:)));

     Tneumann = e([1 2 5], find(~Tneunan)); %Create boundary condition matrices.
     TneumannT = e([1 2 5], find(~TneuTnan));
     Tdirichlet = e([1 2 5], find(~Tdirnan));

     Wneumann = e([1 2 5], find(~Wneunan));
     WneumannT = e([1 2 5], find(~WneuTnan));
     Wdirichlet = e([1 2 5], find(~Wdirnan));

     Uneumann = e([1 2 5], find(~Uneunan));
     UneumannT = e([1 2 5], find(~UneuTnan));
     Udirichlet = e([1 2 5], find(~Udirnan));

     %Create vectors containing the free nodes (non dirichlet nodes)
     TFreeNodes = setdiff(1:pCount,unique(Tdirichlet));
     uFreeNodes = setdiff(1:pCount, unique(Udirichlet));
     wFreeNodes = setdiff(1:pCount, unique(Wdirichlet));

     TFreeCount = size(TFreeNodes,2);
     uFreeCount = size(uFreeNodes,2);
     wFreeCount = size(wFreeNodes,2);
    
     if(iterN == 1 && T0 == (Tref+3))
       uGuess = zeros(uFreeCount,1);
       wGuess = zeros(wFreeCount,1);
       TGuess = Tref*ones(TFreeCount,1);
     else
       uGuess = full(uu(uFreeNodes));
       wGuess = full(wu(wFreeNodes));
       TGuess = full(Tu(TFreeNodes));
     end

     %Create pointers to fields in the results vector.
     uStart = 1;
     uEnd = uFreeCount;
     wStart = uEnd + 1;
     wEnd = wStart + wFreeCount - 1;
     TStart = wEnd + 1;
     TEnd = TStart + TFreeCount - 1;

     uNodes = uStart:uEnd;
     wNodes = wStart:wEnd;
     TNodes = TStart:TEnd;


     %Assemble load vector

     fprintf(1, 'Assembling load vector...');

     for k = 1:tCount
       fw(t(1:3,k),1) = fw(t(1:3,k),1) - g*beta*Tref*area(k)/3;
     end

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
	   alpha*ones(2,1)*TneumannConditions(Tneumann(3,j))/2;
     end


     for j = 1:size(Uneumann,2)
       Gu(Uneumann(1:2,j))=Gu(Uneumann(1:2,j)) + ... 
	   norm(p(:,Uneumann(1,j)) - ... 
		p(:,Uneumann(2,j)))* ...
	   nu*ones(2,1)*UneumannConditions(Uneumann(3,j))/2;
     end

     for j = 1:size(Wneumann,2)
       Gw(Wneumann(1:2,j))=Gw(Wneumann(1:2,j)) + ... 
	   norm(p(:,Wneumann(1,j)) - ... 
		p(:,Wneumann(2,j)))* ...
	   nu*ones(2,1)*WneumannConditions(Wneumann(3,j))/2;
     end
     
     for k = 1:size(TneumannT,2)
       L = norm(p(:,TneumannT(1,k)) - p(:,TneumannT(2,k)));
       
       QT(TneumannT([1 2],k), TneumannT([1 2],k)) = QT(TneumannT([1 ...
		    2],k),TneumannT([1 2],k)) ...
	   - alpha*TneumannTConditions(TneumannT(3,k))*L*[2,1;1,2]/6;
     end
     
     fu = fu + Gu;
     fw = fw + Gw;
     fT = fT + GT;
     
     %Note that we don't have any speed dependance in our Neumann
     %boundary conditions for the velocity vector.
     %Also note that we are setting the divergence of the velocity
     %vector to zero on the border. 
     
     fprintf(1, ' done!\n')


     
     fprintf(1, 'Enforcing dirichlet conditions...')
     uu(:) = 0;
     wu(:) = 0;
     Tu(:) = 0;

     Tu(Tdirichlet(1,:)) = TdirichletConditions(Tdirichlet(3,:));
     Tu(Tdirichlet(2,:)) = TdirichletConditions(Tdirichlet(3,:));
     uu(Udirichlet(1,:)) = UdirichletConditions(Udirichlet(3,:));
     uu(Udirichlet(2,:)) = UdirichletConditions(Udirichlet(3,:));
     wu(Wdirichlet(1,:)) = WdirichletConditions(Wdirichlet(3,:));
     wu(Wdirichlet(2,:)) = WdirichletConditions(Wdirichlet(3,:));


     %Put dirichlet conditions into the solution vectors
     
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
     wu(int(ind)) = 0.5*(WdirichletConditions(Wdirichlet(3,ca(ind))) + ...
			 WdirichletConditions(Wdirichlet(3,cb(ind))));
     
     %Update load vectors with dirichlet conditions by moving them to
     %the RHS.

     fu = fu - LHSu(uu, wu, Tu);    
     fw = fw - LHSw(uu, wu, Tu);
     fT = fT - LHST(uu, wu, Tu);
     fprintf(1, ' done!\n');


     [resu, resw, resT, ~, ~] = newtonraphson(fu, fw, fT, uGuess, ...
					      wGuess, TGuess, uFreeNodes, ...
					      wFreeNodes, ...
					      TFreeNodes,uNodes,wNodes, TNodes, ITERMAX, epsilon);
    

     Tu(TFreeNodes) = resT;
     uu(uFreeNodes) = resu;
     wu(wFreeNodes) = resw;
     edges = find(~(e(5,:)-WallID));
     meantemp = mean(Tu(unique([e(1,edges) e(2,edges)])));
     disp(['MeanTemp: ' num2str(meantemp) ', T0 = ' num2str(T0)]);
     newtemp = meantemp;
  end

  %Calculate the U-value of the air
  %and put the result in a vector

  %Get the blackbody radiation and solar intensity
  %NeumannConditions and NeumannTConditions?
  
  QConvection = kAir*(TneumannConditions(WallID)+meantemp*TneumannTConditions(WallID));
  
  AirH(iterN) = QConvection/(meantemp - Tref);
  disp(['Found H-value: ' num2str(AirH(iterN))])
end

figure(1)
hold off
plot(Range, AirH)
xlabel(VarName)
ylabel('h-value')


ret = [Range', AirH];

figure(2)
hold off
tricontourf(p(1:2,:),t(1:3,:),sqrt(uu.^2+wu.^2))

xlabel('Position (m)')
ylabel('Position (m)')
title('Beloppet av luftens hastighetsvektor (m/s)')

figure(3)
hold off
triquiver(p(1:2,:), t(1:3,:), uu, wu,20, [0 20], [0 22])
hold on
plot([0,20,20,-6,0,0], [0,0,22,22,18,0], 'k')
hold off
xlim([-6 20])
ylim([0 22])
xlabel('Position (m)')
ylabel('Position (m)')
title('Vindriktning')
%figure(3)
%hold off
%sCount = 3;
%tristreamline(p(1:2,:),t(1:3,:), uu, wu, [5], [16]);
%triquiver(p(1:2,:), t(1:3,:), uu, wu,20, [3 17], [2 20])
 %(22*[0:(sCount-1)]/(sCount-1)+1)')
figure(4)
hold off
tricontourf(p(1:2,:),t(1:3,:),Tu-kelvin)

xlabel('Position (m)')
ylabel('Position (m)')
title('Temperatur (C)')