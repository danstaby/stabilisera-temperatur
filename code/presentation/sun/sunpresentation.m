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
hline=refline(0,0);
set(hline,'Color','k')
xlim([0,24])
set(gca,'XTick',[0:4:24],'XTickLabel',[0:4:24])
xlabel('Klockslag, UTC+1')
ylabel('Flöden genom fönster, Wm^{-2}')

gtext('Solens intensitet')
gtext('Solinstrålning')
gtext('Värmeledning och konvektion')
gtext('Långvågsstrålning')