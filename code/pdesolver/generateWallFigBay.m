function ret = generateWallFigBay(fileName)

load(fileName);
timeOffset = 28*24;

ind = [28*24*3600/500-0.4:29*24*3600/500-0.2];
time = (data1(ind,1)-data1(ind(1),1))/3600;

hold off
fig1 = figure(1);
plot(time, data1(ind,2))
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
title('Kyleffekt utskjutande fasaddelar', 'fontsize', 12)
xlim([0 24])
set(fig1, 'Position', [300 700 300 200])
disp('1: Kyleffekt en solig aprildag')

fig2 = figure(2);
plot(time,data2(ind,2))
title('Kyleffekt utskjutande fasaddelar', 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])

set(fig2, 'Position', [800 700 300 200])
disp('2: Kyleffekt en molning aprildag')

fig3 = figure(3);
plot(time,data3(ind,2))
title('Kyleffekt utskjutande fasaddelar', 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])
disp('3: Kyleffekt en molnig decemberdag')
set(fig3, 'Position', [1130 700 300 200])
