% Kod för att testa följande funktioner:
% sunposition(long, lat, y, m, d, hr)
% angletheta(alpha, beta, gamma)
% gvalue(g0, p, q, theta)
% effekt(I, month, day, hour)

clc
close all

long = 11.979435;
lat = 57.691522;

year = 2012;
month = 4;
day = 1;

%% Test av sunposition() och angletheta()

el = [];
az = [];
theta = [];
n = 0;
xax=[];
for hr = 0:23
    for min = 0:59
        UTC = hr + min/60;
        n = n + 1;
        [el(n), az(n)] = sunposition(long, lat, year, month, day, UTC);
        theta(n) = angletheta(el(n), az(n), 180);
        xax(n)=n;
    end
end

xax =xax./60;
plot(xax, theta)
hold on
plot(xax, el, 'r', xax, az, 'g')
legend('Relative to horisontol southern direction','Elevation relative horizon','Azimuth relative eastern direction')
xlabel('Time, UTC')
ylabel('Angle, degrees')

%% Test av effekt()

% Skapa vektor med intensitet I
C = [-0.025 0.6 -2.7];
I = [];
m = 0;
for hour = 0:23
    for min = 0:59
        m = m + 1;
        x = hour + min/60;
        if x > 6 & x < 18
            I(m) = C(1)*x^2 + C(2)*x + C(3);
        else
            I(m) = 0;
        end
    end
end
%figure(2)
%plot(1:m, I)

% Test
eff = [];
o = 0;
xax2=[];
I0=200;
for hour = 0:23
    for min = 0:59
        o = o + 1;
        t = hour + min/60;
        eff(o) = effekt(I0, month, day, t);
        xax2(o)=o;
    end
end

xax2=xax2./60;
figure(2)
plot(xax2, eff)
ylabel('Effekt, W')
xlabel('Tid, UTC')