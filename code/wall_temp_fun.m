% IN:    -utomhustemperatur och solinstr�lning
% OUT: -V�ggens temperatur
%
% Kan �ven ber�kna energifl�det genom v�ggen.
% samt hur l�ng tid det tar f�r v�ggen att komma upp i temp.


function [output] = wall_temp_fun(T_u, I)

% KONSTANTER

L=0.5; %m, v�ggens tjocklek
lambda=0.593; %W/mK, tegels v�rmeledningsf�rm�ga
U1=lambda/L;

sigma=5.67*10^(-8); %W/m2K4
kp=50; % konvektionsparametern
abs=0.8; % V�ggens absorptionsparameter f�r solstr�lning

T_i=20+273; % K, inomhustemperatur
T_u=T_u;
T_a=T_u-30; %K atmosf�rens temperatur

Cv=0.8; % W/kgK
rho=1600; %kg/m3
m=L*rho; % massa kg/m2
C=Cv*m; % Ws/Km2

% I=600; %W/m2 Instr�lad solintensitet

% V�GGTEMPTERTUREN vid konstant ute- och innetemperatur
a=-sigma/C;
b=(-kp+U1)/C;
c=((0.2.*T_u.^4+0.8.*T_a.^4).*sigma+abs.*I+kp.*T_u-T_i*U1)/C;

T_w=(0:0.01:30)+273; % i C
y_w=a.*T_w.^4+b.*T_w+c; %y=dT/dt, y=0

% s�k STEADYSTATE (nollst�lle t. y_w)
y_w=y_w.^2;
miny=min(y_w);
index_w=find(y_w<=miny);
zero_Tw=T_w(index_w);


% V�GGTEMPTERTUREN isolerad v�gg
L_i=0.1; %m, isoleringens tjocklek
lambda_i=0.037; %W/mK, mineralulls v�rmeledningsf�rm�ga
U_i=lambda_i/L_i;
U2=1/(1/U1+1/U_i);

a=-sigma;
b=(-kp+U2);
c=(0.2.*T_u.^4+0.8.*T_a.^4).*sigma+abs*I+kp.*T_u-T_i*U2;

T_w=(0:0.01:30)+273; % i C
y_w=a.*T_w.^4+b.*T_w+c; %y=dT/dt, y=0

% s�k STEADYSTATE (nollst�lle t. y_w)
y_w=y_w.^2;
miny=min(y_w);
index_w=find(y_w<=miny);
zero_Tw2=T_w(index_w);



% Fl�det ut ur en oisolerad s�derv�gg per kvadr.meter
F1=U1*(T_i - zero_Tw);

% Fl�det ut ur en isolerad typ norrv�gg per kvadr.meter
F2=U2*(T_i - zero_Tw2);

% Fl�det ut ur en genomsnittlig kvadratmeter v�gg, med oisolerad s�derv�gg
U3=0.6086;
F3=U3*(T_i - zero_Tw);

% Fl�det ut ur en genomsnittlig kvadratmeter v�gg, med isolerad s�derv�gg
U4=0.4680;
F4=U4*(T_i - zero_Tw2);

% Fl�det ut ur en genomsnittlig kvadratmeter v�gg, med isolerad s�der- och v�sterv�gg
U5=0.4112;
F5=U5*(T_i - zero_Tw2);


% S� l�ng tid tar det f�r v�ggens yta att komma upp i temperatur.
% t_utan=int_solution(T_u, zero_Tw);
% t_mediso=int_solution(T_u, zero_Tw2);

output=[zero_Tw];




