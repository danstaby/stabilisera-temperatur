function Q = effekt(I, month, day, hour)
% hour in UTC decimal format
% 0 < I < 1400 W/m^2, typ
A = 1; % Sammanlagd area av sydsidans fönster
g0 = 0.61; % Från ASHRAE
p = 3; % Antal glasskivor i fönstret
q = 4; % Beroende av beläggningar o.d., inga belägningar => 4

long = 11.98;
lat = 57.7;
year = 2012;

[alpha, beta] = sunposition(long, lat, year, month, day, hour);

theta = angletheta(alpha, beta, 180);

g = gvalue(g0, p, q, theta);

Q = g*A*I*cosd(theta);