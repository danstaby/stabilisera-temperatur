function Q = effekt(I, month, day, hour)
% hour in UTC decimal format

% 0 < I < 1400 W/m^2, typ

A = 1.5; % Fonstrets area
g0 = 0.55;
p = 3; % Antal glasskivor i fonstret
q = 4; % Beroende av belaggningar o.d., inga belagningar => 4

long = 11.98;
lat = 57.7;
year = 2012;

[alpha, beta] = sunposition(long, lat, year, month, day, hour);

theta = angletheta(alpha, beta, 180)

g = gvalue(g0, p, q, theta);

Q = g*A*I;