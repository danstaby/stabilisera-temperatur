function ret = generateWallFigures(fileName)

load(fileName);
timeOffset = 28*24;

ind = [28*24*3600/500-0.4:29*24*3600/500-0.2];
time = (aprsun1(ind,1)-aprsun1(ind(1),1))/3600;

makePlot(time, aprsun1, ind, 'tegel', 1);
makePlot(time, aprsun2, ind, 'tegel och isolering', 2);
makePlot(time, aprsun3, ind, 'utskjutande fasaddelar', 3);
makePlot(time, aprnosun1, ind, 'tegel', 4);
makePlot(time, aprnosun2, ind, 'tegel och isolering', 5);
makePlot(time, aprnosun3, ind, 'utskjutande fasaddelar', 6);
makePlot(time, decsun1, ind, 'tegel', 7);
makePlot(time, decsun2, ind, 'tegel och isolering', 8);
makePlot(time, decsun3, ind, 'utskjutande fasaddelar', 9);
makePlot(time, decnosun1, ind, 'tegel', 10);
makePlot(time, decnosun2, ind, 'tegel och isolering', 11);
makePlot(time, decnosun3, ind, 'utskjutande fasaddelar', 12);

function makePlot(time,data, ind, Type, ID)

fig = figure(ID);
hold off
plot(time,data(ind,2))
title(['Kyleffekt, ' Type], 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])
set(fig, 'Position', [303*(mod(ID-1,4)) 700-(floor(ID/4))*203 300 200])