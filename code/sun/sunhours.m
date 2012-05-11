function ret = sunhours()

strMonths = {'jan', 'feb', 'mar', 'apr', 'maj', 'jun', 'jul', 'aug', ...
	    'sep', 'okt', 'nov', 'dec'};
Warm = [9, 12+4];

avg = (Primitive(Warm(2))-Primitive(Warm(1)))/...
      (Warm(2)-Warm(1));

months = Warm(1):0.05:Warm(2);
h = figure(1);
plot(months, TotalHours(months),'r')
hold on
plot(months, avg*ones(1,size(months,2)), '--', 'Color', 'k')
ylabel('Soltimmar per dygn')
xlabel('Manad')
set(gca, 'Xtick', Warm(1):Warm(2), ...
	 'XTickLabel', strMonths(mod((Warm(1):Warm(2)),12)+1))
hold off
legend('Soltimmar', 'Medel', 'Location', 'best')
set(h, 'Position', [100,100,400,200])
xlim(Warm)
ret = avg;

function val = TotalHours(t)
val = 12-6*cos(pi*(t+1/3)/6);

function val = Primitive(t)
val = 12*t-36*sin(pi*(t+1/3)/6)/pi;