%clc
%clear all


% Beräknar solens intensitet över dygnet.

long=12.0; % grader
lat=57.7; % grader
y=2011;
r=6371.2; % m, jordens radie
I0=1370; % W/m2
mu=4.6*10^(-5); % /m


% DECEMBER
m=12;
d=31;
for i=1:240
    hr=0.1*i;
    [alpha beta]=sunposition(long, lat, y, m, d, hr);
    xdec(i)=r*cosd(alpha+90)+sqrt((r*cosd(alpha+90))^2-r^2+(r+15000)^2);
    Idec(i,1)=I0*exp(-mu*xdec(i));
    Idec(i,2)=hr;
    if alpha < 0
        Idec(i,1) = 0;
    end
end

Idec(:,2)=Idec(:,2)+2;
for i=1:length(Idec)
    if Idec(i,2)>24
       Idec(i,2)=Idec(i,2)-24;
    end
end

