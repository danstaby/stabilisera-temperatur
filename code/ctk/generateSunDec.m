function ret = generateSunDec()

load('sunpower.mat');
load('rawwalldata.mat');
load('roofheat.mat');




timeOffset = 28*24;

ind = [28*24*3600/500-0.4:29*24*3600/500-0.2];
time = (decnosun1(ind,1)-decnosun1(ind(1),1))/3600;


ANorth = 290;
AOther = 151+61;
AWindow = 109;
ARoof = 257/2;
ABay = 47;

Energy = [];

%Roof contribution
% - sunnorth
% - sunsouth
% - nosun
figure(2)
Energy = -(sunnorth(ind,2)-nosun(ind,2))*ARoof ...
       - (sunsouth(ind,2)-nosun(ind,2))*ARoof;

plot(time, -(sunnorth(ind,2)-nosun(ind,2))*ARoof ...
       - (sunsouth(ind,2)-nosun(ind,2))*ARoof)

hold on

%Bay wall contribution
Energy = Energy - (decsun3(ind,2)-decnosun3(ind,2))*ABay;

plot(time, -(decsun3(ind,2)-decnosun3(ind,2))*ABay, 'r')

%Other wall contribution

Energy = Energy - (decsun1(ind,2) - decnosun1(ind,2))*AOther;
plot(time, -(decsun1(ind,2) - decnosun1(ind,2))*AOther, 'g')
%North wall contribution

Energy = Energy - (decsun2(ind,2) - decnosun2(ind,2))*ANorth;
plot(time, -(decsun2(ind,2) - decnosun2(ind,2))*ANorth, 'k')
%Window contribution

legend('Roof', 'Bay', 'Other', 'North')
%Energy = Energy + AWindow*sunpowerdec(:,2);
plot(time, AWindow*sunpowerdec(:,2), 'm')
hold off
legend('Roof', 'Bay', 'Other', 'North', 'Window')
figure(1)
plot(time, Energy)
hold on
xlim([0 24]);
plot(time, mean(Energy)*ones(size(time,1),2), '--', 'Color', ...
	 'k')
legend('Momentan', 'Medel', 'location', 'best')
xlabel('Tid (h)')
ylabel('Effekt (W)')
mean(Energy)
hold off

figure(3)
plot(time, sunpowerdec(:,2))
xlim([0 24])

%{
hold off
fig1 = figure(1);
plot(time, data1(ind,2))
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
title('Kyleffekt med enbart tegel', 'fontsize', 12)
xlim([0 24])
set(fig1, 'Position', [300 700 300 200])


fig2 = figure(2);
plot(time,data2(ind,2))
title('Kyleffekt med tegel och isolering', 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])

set(fig2, 'Position', [800 700 300 200])
%}