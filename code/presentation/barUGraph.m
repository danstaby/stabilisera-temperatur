function ret = barUGraph()
clc
str=sprintf('Syd- och\n\r v�stv�gg');
U = [1.186,0.279,0.393, 0.171, 1];
strNames = {str, ...
	    'Norrv�gg',...
	    'Burspr�k',...
	    'Tak',...
	    'F�nster'};
cvec = [2,1,6,4,7];


cm = colormap(hsv(7));

% 1 - norr
% 2 - SV
% 6 - bay
% 4 - Tak
% 5 - grund
% 3 - konstnta
% 7 - F�nster


figure(1)
hold off
data = NaN(5,5);
for n = 1:5
  data(n,n) = U(n);
end

bh = bar(data, 'stacked');


for n = 1:max(size(bh))
  set(bh(n), 'FaceColor', cm(cvec(n),:));
end

ylabel('U-v�rde (Wm^{-2}s^{-1})', 'Fontsize', 14);
title('Klimatskalets U-v�rden', 'Fontsize', 22);
ax = axis; % Current axis limits
axis(axis); % Set the axis limit modes (e.g. XLimMode) to manual
Yl = ax(3:4); % Y-axis limits

% Place the text labels
t = text(1:5,zeros(1,5),strNames);
set(t,'HorizontalAlignment','left','VerticalAlignment','top', ...
'Rotation',-45, 'Units', 'normalized', 'Clipping', 'off', 'Fontsize', 14);

% Remove the default labels
set(gca,'XTickLabel','')

TickExt = zeros(5,4);
for ii = 1:5
    TickExt(ii,:) = get(t(ii),'Extent');
end
Ha = gca;
Hfig = get(Ha, 'Parent');
set(Hfig, 'Units', 'normalized');
set(Ha, 'Units', 'normalized', 'ActivePositionProperty', 'OuterPosition');
Axpos = get(Ha, 'OuterPosition');
%Move axis object

TickExt;

offs = min(TickExt(:,2))-max(TickExt(:,2))+0.05;%-TickExt(:,4))
Axpos(2) = Axpos(2)-offs;
Axpos(4) = Axpos(4)+offs;

set(Ha,'OuterPosition',Axpos);

set(Hfig, 'Units', 'normalized');
set(Ha, 'Units', 'normalized');
set(t, 'Units', 'normalized');

U = [1.186,0.279,0.393, 0.171, 1];
% Add values for bars
p1 = text(0.8,1.25,'1,19');
p2 = text(1.78,0.35,'0,28'); 
p3 = text(2.8,0.45,'0,39');
p4 = text(3.8,0.23,'0,17');
p5 = text(4.8,1.06,'1,00'); 
set(p1, 'Fontsize', 14);
set(p2, 'Fontsize', 14);
set(p3, 'Fontsize', 14);
set(p4, 'Fontsize', 14);
set(p5, 'Fontsize', 14);




