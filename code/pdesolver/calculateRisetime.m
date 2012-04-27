function ret = calculateRisetime(fileName)

load(fileName);

%step1 is the wall without insulation and
%step 2 is with extra insulation



%Norm data

norm1 = normData(step1);
norm2 = normData(step2);

t10 = [0,0];
t90 = [0,0];

%Find 10% point
t10(1) = findTime(norm1, 0.1);
t10(2) = findTime(norm2, 0.1);

%Find 90% point

t90(1) = findTime(norm1, 0.9);
t90(2) = findTime(norm2, 0.9);

risetime = (t10-t90)/3600;

%Display data

disp(['Rise time of brick wall: ' num2str(risetime(1))])
disp(['Rise time of insulated wall: ' num2str(risetime(2))])

%Plot data
makePlot(step1, t10(1), t90(1), 1, 'Enbart tegel');
makePlot(step2, t10(2), t90(2), 2, 'Tegel och mineralull');



function ret = makePlot(data, t10, t90, figId, strTitle)
tOffset = 3600;
figure(figId)
hold off
plot(data(:,1)/tOffset, data(:,2))
hold on
yl = ylim;

plot([t90, t90]/tOffset, [yl(1)-1, yl(2)+1],'--','Color','k')
plot([t10, t10]/tOffset, [yl(1)-1 yl(2)+1],'--', 'Color', 'k')
hold off
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
title(strTitle, 'fontsize', 12)
xlim([0 t10/tOffset+24*3600/tOffset])
ylim([min(data(:,2)), max(data(:,2))])
function normed = normData(data)

Qmax = max(data(:,2));
Qmin = min(data(:,2));
normed = data;
normed(:,2) = (data(:,2)-Qmin)/(Qmax-Qmin);



function time = findTime(data, percentage)
%findTime(data, percentage)
%
%This function only works for decreasing data-sets


%Find points
ind = find(data(:,2) < percentage);
p1 = ind(1)-1;
p2 = ind(1);

%Make spline interpolation and return data
time = data(p1,1) + (percentage-data(p1,2))*(data(p2,1)-data(p1,1))/(data(p2,2)-data(p1,2))

