function ret = sunfigures(month, day, date)
close all

if month>3 month<10
    sommartid=1;
else
    sommartid=0;
end

year = 2011;

long = 11.979435;
lat = 57.691522;

r = 6.731e6;
I0 = 1370;
mu = 4.6e-5;

el = [];
az = [];
theta = [];
x = [];
I = [];
eff = [];
n = 0;
xax=[];

for hr = 0:23
    for min = 0:59
        UTC = hr + min/60;
        n = n + 1;
        [el(n), az(n)] = sunposition(long, lat, year, month, day, UTC);
        theta(n) = angletheta(el(n), az(n), 180);
        if el(n) > 0
            x(n) = r*cosd(90+el(n))+sqrt((r*cosd(90+el(n)))^2+((30e3)*r+(15e3)^2));
            I(n) = I0*exp(-mu*x(n));
        else
            x(n) = 0;
            I(n) = 0;
        end
        eff(n) = effekt(I(n), month, day, UTC);
        xax(n)=n;
    end
end

figure(1)
xax =xax./60+1+sommartid;
plot(xax, theta)
hold on
plot(xax, el, 'r', xax, az, 'g')
xlabel('Tid, UTC+2 [timmar]')
ylabel('Vinkel [grader]')
axis([1+sommartid 22+1+sommartid -100 500])
legend('Relativt horisontella södra riktningen','Höjd över horisonten','Azimuthala relativt öster, medsols positivt')
%legend(strcat('Relativt horisontella södra riktningen, ', num2str(date)), ...
%    strcat('Höjd över horisonten, ', num2str(date)), ...
%    strcat('Azimuthala relativt öster, medsols positivt, ', num2str(date)))

figure(2)
plot(xax, I)
hold on
plot(xax, eff, 'r')
xlabel('Time, UTC [h]')
ylabel('Intensity [W/m^2]')
axis([0 24 0 1370])
legend(strcat('Suns intensity, ', num2str(date)), ...
    strcat('Flux through window, ', num2str(date)))

hold on

figure(3)
plot(xax, x)

efftotal = (1/60).*trapz(eff);
I15 = I0*exp(-mu*15e3)
ret = efftotal;