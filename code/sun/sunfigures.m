function [timefinal, intensity] = sunfigures(month, day, date)
% Month and day integers, date is string, output is a 2 column matrix
close all

monthvector = [31,28,31,30,31,30,31,31,30,31,30,31];
currentmonths = zeros(1,12);
for i = 1:(month-1)
      currentmonths(i)=1;
end

daynumber = dot(monthvector,currentmonths)

if month>3 month<10
    sommartid=1;
else
    sommartid=0;
end

year = 2011;

long = 11.979435;
lat = 57.691522;

I0 = 1370;
tau = [0.142,0.144,0.156,0.180,0.196,0.205,0.207,0.201,0.177,0.160,0.149, 0.142];
daynumber = 360*(daynumber-1)/365;
Ispace = I0*(1.00011 + 0.034221*cosd(daynumber) + 0.00128*sind(daynumber) + ...
	  0.000719*cosd(2*daynumber) + 0.000077*sind(daynumber));
el = [];
az = [];
theta = [];
I = [];
eff = [];
n = 0;

for time=-(1+sommartid)*3600:500:((23-sommartid)*3600+100);
    UTC = time/3600;
    n = n + 1;
    if time>0
        [el(n), az(n)] = sunposition(long, lat, year, month, day, UTC);
        theta(n) = angletheta(el(n), az(n), 180);
    else
        el(n) = 0;
        az(n) = 0;
        theta(n) = 0;
    end

    if el(n) > 0
        I(n) = Ispace*exp(-tau(month)/cosd(90-el(n)));
    else
        I(n)=0;
    end
    eff(n) = effekt(I(n), month, day, UTC);
    x(n)=UTC;
    
end

%x = x/3600;
x = x + 1 + sommartid;

figure(1)
plot(x, theta)
hold on
plot(x, el, 'r', x, az, 'g')
xlabel('Tid, UTC+2 [timmar]')
ylabel('Vinkel [grader]')
axis([1+sommartid 22+1+sommartid -100 500])
legend('Relativt horisontella södra riktningen','Höjd över horisonten','Azimuthala relativt öster, medsols positivt')

figure(2)

plot(x, I, x, eff, 'r')
xlabel('Tid, UTC [h]')
ylabel('Intensitet [Wm^{-2}]')
axis([0 24 0 1370])
legend(['Solens intensitet, ' date], ...
       ['Intensitet genom glaset, ' date])

efftotal = (1/60).*trapz(eff);

timefinal = x;
intensity = I;
