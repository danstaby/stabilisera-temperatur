function ret = calculatePhase(fileName)


load(fileName);
timeOffset = 28*24;

norm1 = normData(data1);
norm2 = normData(data2);

[ho1, peakout1] = findpeaks(norm1(:,2),'NPEAKS', 1)
[hi1, peakin1] = findpeaks(norm1(:,3), 'NPEAKS', 1)
(norm1(peakout1,1)-norm1(peakin1,1))/3600

figure(1)
plot(norm1(:,1)/(3600)-timeOffset,norm1(:,2))
hold on
plot(norm1(:,1)/(3600)-timeOffset,norm1(:,3), 'r')

plot(norm1(peakout1,1)/3600-timeOffset, ho1, 'o')
plot(norm1(peakin1,1)/3600-timeOffset, hi1, 'o')
title('Normerad energi utan isolering')
xlabel('Tid i timmar')
ylabel('Normerad energi (Ingen enhet)')
legend('Energi in', 'Energi ut')
xlim([0 24])
hold off

figure(2)
plot(norm2(:,1)/(3600)-timeOffset,norm2(:,2))
hold on
plot(norm2(:,1)/(3600)-timeOffset,norm2(:,3), 'r')
title('Normerad energi med isolering')
xlabel('Tid i timmar')
ylabel('Normerad energy')
legend('Energi in', 'Energi ut')
hold off
xlim([0 24])

function ret = normData(rawData)

int = [(28*24*3600/500-0.4):(29*24*3600/500-0.2)];

ret = rawData(int,:);
ret(:,2) = (ret(:,2)-mean(ret(:,2)))/(max(ret(:,2)-mean(ret(:,2))));
ret(:,3) = (ret(:,3)-mean(ret(:,3)))/(max(ret(:,3)-mean(ret(:,3))));