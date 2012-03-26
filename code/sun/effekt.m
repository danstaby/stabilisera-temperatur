function Q = effekt(I, hour)

% 0 < I < 1300 W/m^2

A = 1.5; % Fönstrets area
g0 = 0.55;
p = 3; % Antal glasskivor i fönstret
q = 4; % Beroende av beläggningar o.d.

long = 11.98;
lat = 57.7;
year = 2012;
month = 6;
day = 25;
% hour in UTC decimal format

[alpha, beta] = sunposition(long, lat, year, month, day, hour);

theta = angletheta(alpha, beta, 90)

g = gvalue(g0, p, q, theta);

Q = g*A*I;