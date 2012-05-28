clc
close all
load('windowsimulation.mat')
load('sunintensitydec.mat')
load('windowsconduction.mat')
load('windowslongwave.mat')

cm=colormap(hsv(7));
windowcolor=hsv(7);

hold on
plot(windowssundec(:,1), windowssundec(:,2), 'Color', [cm(7,:)])
plot(intensitydec(:,1), intensitydec(:,2), '--', 'Color', [cm(7,:)])
plot(conductionwindows(:,1), conductionwindows(:,2), 'Color', [cm(7,:)])
plot(longwavewindows(:,1),longwavewindows(:,2), 'Color', [cm(7,:)])
%hline=refline(0,0);

plot(windowssundec(:,1),zeros(1,174),'--k')
%set(hline,'Color','k')
xlim([0,24])
set(gca,'XTick',[0:4:24],'XTickLabel',[0:4:24])
%title('Fl�den genom f�nster', 'Fontsize', 22)
xlabel('Tid, timmar', 'Fontsize', 14)
ylabel('Energifl�den, Wm ^{-2}', 'Fontsize', 14)

gtext('Solens intensitet', 'Fontsize', 14)
gtext('Solinstr�lning', 'Fontsize', 14)
gtext('V�rmeledning och konvektion', 'Fontsize', 14)
gtext('L�ngv�gsstr�lning', 'Fontsize', 14)