function ret = generateWallFigNew(fileName)

load(fileName);
timeOffset = 28*24;

ind = [28*24*3600/500-0.4:29*24*3600/500-0.2];
time = (decsun1(ind,1)-decsun1(ind(1),1))/3600;

hold off
fig1 = figure(1);
plot(time, decsun1(ind,2))
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
title('Kyleffekt med enbart tegel', 'fontsize', 12)
xlim([0 24])
set(fig1, 'Position', [300 700 300 200])


fig2 = figure(2);
plot(time,decsun2(ind,2))
title('Kyleffekt med tegel och isolering', 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])

set(fig2, 'Position', [800 700 300 200])

fig3 = figure(3);
plot(time,decsun3(ind,2))
title('Kyleffekt utskjutande fasaddelar', 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])

set(fig3, 'Position', [1150 700 300 200])

