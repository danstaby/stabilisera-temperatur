function [alpha, beta] = sunposition(long, lat, y, m, d, hr)
%{
long positive east of Greenwich, in degrees
lat positive north of equator, in degrees
y = year
m = month
d = day
hr = hour in UTC decimal form, including minute and seconds
alpha = elevation
beta = azimuth
%}

%{
For Walleriusgatan:
long = 11.979435;
lat = 57.691522;
%}

% Degrees to radians
rad = pi/180;
% Earth mean radius, km
emr = 6371.01;
% Astronomical unit, km
au = 149597890;

% Calculate julian day
% Integer division for everything except hr
jd = floor((1461*(y+4800+floor((m-14)/12)))/4) ...
    + floor((367*(m-2-12*(floor((m-14)/12))))/12) ...
    - floor((3*((y+4900+floor((m-14)/12))/100))/4) ...
    + d - 32075 - 0.5 + hr/24;

% Difference betwween current day and noon 1 jan 2000 (UTC)
n = jd - 2451545;

% Ecliptic coordinates, all angles in radians
Omega = 2.1429 - 1.0394594e-3*n;
% Mean longitude
L = 4.895063 + 0.017202791698*n;
% Mean anomaly
g = 6.24006 + 0.0172019699*n;
% Ecliptic longitude
l = L + 0.03341607*sin(g) + 3.4894e-4*sin(2*g) ...
    - 1.134e-4 - 2.03e-5*sin(Omega);
% Obliquity of ecliptic
ep = 0.4090928 - 6.214e-9*n + 3.96e-5*cos(Omega);

% From ecliptic to celestial coodinates
% Right ascension (must be larger than 0)
ra = atan2((cos(ep)*sin(l)), cos(l));
if ra < 0
    ra = ra + 2*pi;
end
% Declination
delta = asin(sin(ep)*sin(l));

% From celestial to horizontol coordinates
% Greenwich mean sidereal time
gmst = 6.6974243242 + 0.0657098283*n + hr;
% Local mean siderial time
lmst = (gmst*15 + long)*rad;
% Hour angle
omega = lmst - ra;
% Latitude in radians
latrad = lat*rad;
% Zenith angle
thetaz = acos(cos(latrad)*cos(omega)*cos(delta) ...
    + sin(delta)*sin(latrad));

% Azimuth angle, positive east of south
Y = -sin(omega);
X = tan(delta)*cos(latrad) - sin(latrad)*cos(omega);
gamma = atan2(Y, X);
if gamma < 0
    gamma = gamma + 2*pi;
end

% Parallax
parallax = emr*sin(thetaz)/au;
% Parallax correction to zenith angle
thetaz = thetaz + parallax;

% Elevation from horizon
theta = pi/2 - thetaz;

% Output in degrees, alpha = elevation (from vertical), beta = azimuth
alpha = theta/rad;
beta = gamma/rad;
