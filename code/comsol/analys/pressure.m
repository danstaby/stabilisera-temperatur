function ret = pressure(file)

TOLERANCE = 1e-6;
MAX_ITER = 40;
pref = 1.01e5;
guess = pref;

data = importdata(file);

wind = data(:,1);
house1 = data(:,2:5);
house2 = data(:,6:9);

walls = [18, 0, 0, 18]';

%Calculate pressure inside the house

x1 = guess*ones(size(wind,1),1);
x2 = x1;

fillvector = ones(1,size(house1,2));

res = TOLERANCE + 1;
n = 0;
while max(abs(res)) > TOLERANCE && n < MAX_ITER
  %Calculate value

  % Sum A_n dp_n^0.6 = 0 
  
  res1 = ((house1-x1*fillvector))*walls;
  res2 = ((house2-x2*fillvector))*walls;

  %Calculate jacobian
  
  %d/dp = Sum -0.6*A_n dp_n^-0.4
  
  %Jac1 = -0.6*((house1-x1*fillvector)).^(-0.4)*walls;
  %Jac2 = -0.6*((house2-x2*fillvector)).^(-0.4)*walls;
  Jac1 = - sum(walls);
  Jac2 = - sum(walls);

  dx1 =  res1./Jac1;
  dx2 =  res2./Jac2;

  %Update value
  x1 = x1 - dx1;
  x2 = x2 - dx2;

  res = [res1; res2];
  n = n + 1;
  disp(['Maximum residue: ' num2str(max(res)) ' Max step: ' ...
       num2str(max([dx1;dx2]))])
end


%Display data
figure(1)
plot(wind, abs(pref - x1))
hold on
plot(wind, abs(pref - x2), 'r')
plot(wind, 0.25*1.2*wind.^2, 'k')
hold off
xlabel('Wind speed (m/s)')
ylabel('Pressure difference (Pa)')
legend('House 1', 'House 2', 'Theoretical')

figure(2)
plot(wind, house1(:,1)-x1)
hold on
plot(wind, house2(:,4)-x2, 'r')
plot(wind, (0.25*1.2*wind.^2).^0.6, 'k')
hold off
xlabel('Wind speed (m/s)')
ylabel('Leakage driving number (Pa)isch')
legend('House 1', 'House 2', 'Theoretical')


house1(:,1)-x1
house2(:,4)-x2

figure(3)
plot(wind, house1(:,1))
hold on
plot(wind, x1 ,'r')
plot(wind,house2(:,4) ,'m')
plot(wind,x2 ,'k')
hold off
xlabel('Wind (m/s)')
ylabel('Pressure (Pa)')
legend('Boundary Pressure House 1', 'Pressure inside house 1', ...
       'Boundary Pressure House 2', 'Pressure inside house 2')

Area = 198;
Cl = 1.2e-3/(50)^0.6;
Cpair = 1005;
rhoair = 1.293;

lk = Cpair*rhoair*Cl;
fig = [];

fig(1) = figure(4);
plot(wind, lk*(house1(:,1)-x1))
hold on
plot(wind, lk*(house1(:,1)-x1).^0.6,'r')

hold off
xlabel('Vindhastighet ortogonalt mot huset (ms^{-1})')
ylabel('Kyleffekt (Wm^{-2}K^{-1})')
title('Hus i lovart')
legend('Darcys lag','Exponent 0,6', 'location', 'best')
xlim([0 24])

fig(2) = figure(5);
plot(wind, lk*(house2(:,4)-x2))
hold on
plot(wind, lk*(house2(:,4)-x2).^0.6,'r')
hold off
xlabel('Vindhastighet ortogonalt mot huset (ms^{-1})')
ylabel('Kyleffekt (Wm^{-2}K^{-1})')
legend('Darcys lag','Exponent 0,6','location','best')
title('Vindskyddat hus')
xlim([0 24])

fig(3) = figure(6);
plot(wind, lk*(0.25*rhoair*wind.^2))
hold on
plot(wind, lk*(0.25*rhoair*wind.^2).^0.6, 'r')
hold off
xlabel('Vindhastighet ortogonalt mot huset (ms^{-1})')
ylabel('Kyleffekt (Wm^{-2}K^{-1})')
legend('Darcys lag','Exponent 0,6', 'location', 'best')
title('Teoretisk approximation')
xlim([0 24])

for n = 1:3
  set(fig(n), 'position', [20+310*(n-1), 400, 300, 200]);
  
end 