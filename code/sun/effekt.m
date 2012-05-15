function q = effekt(I, month, day, hour)
% hour in UTC decimal format
% I = intensity
g0 = 0.61; % From ASHRAE
p = 3; % Number of glazings
q = 4; % Depending on coatings etc. No coatings => 4
gamma = 180-32; % Angle of windows normal relative north

long = 11.979435;
lat = 57.691522;
year = 2012;

[alpha, beta] = sunposition(long, lat, year, month, day, hour);

theta = angletheta(alpha, beta, gamma);

g = gvalue(g0, p, q, theta);

q = g*I*cosd(theta);
