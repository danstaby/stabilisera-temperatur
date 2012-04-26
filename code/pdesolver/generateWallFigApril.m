function ret = generateWallFigApril(fileName)

load(fileName);
timeOffset = 28*24;

ind = [28*24*3600/500-0.4:29*24*3600/500-0.2];
time = (data1(ind,1)-data1(ind(1),1))/3600;

hold off
figure(1)
plot(time, data1(ind,2))
xlabel('Tid i timmar')
ylabel('Kyleffekt (W/m^2)')
title('Kyleffekt med enbart tegel')
xlim([0 24])

figure(2)
plot(time,data2(ind,2))
title('Kyleffekt med tegel och isolering')
xlabel('Tid i timmar')
ylabel('Kyleffekt (W/m^2)')
xlim([0 24])