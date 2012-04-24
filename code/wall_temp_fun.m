% IN:    -utomhustemperatur och solinstrålning
% OUT: -Väggens temperatur
%
% Kan även beräkna energiflödet genom väggen.
% samt hur lång tid det tar för väggen att komma upp i temp.


function [output] = wall_temp_fun(T_u, I)

% KONSTANTER

L=0.5; %m, väggens tjocklek
lambda=0.593; %W/mK, tegels värmeledningsförmåga
U1=lambda/L;

sigma=5.67*10^(-8); %W/m2K4
kp=50; % konvektionsparametern
abs=0.8; % Väggens absorptionsparameter för solstrålning

T_i=20+273; % K, inomhustemperatur
T_u=T_u;
T_a=T_u-30; %K atmosfärens temperatur

Cv=0.8; % W/kgK
rho=1600; %kg/m3
m=L*rho; % massa kg/m2
C=Cv*m; % Ws/Km2

% I=600; %W/m2 Instrålad solintensitet

% VÄGGTEMPTERTUREN vid konstant ute- och innetemperatur
a=-sigma/C;
b=(-kp+U1)/C;
c=((0.2.*T_u.^4+0.8.*T_a.^4).*sigma+abs.*I+kp.*T_u-T_i*U1)/C;

T_w=(0:0.01:30)+273; % i C
y_w=a.*T_w.^4+b.*T_w+c; %y=dT/dt, y=0

% sök STEADYSTATE (nollställe t. y_w)
y_w=y_w.^2;
miny=min(y_w);
index_w=find(y_w<=miny);
zero_Tw=T_w(index_w);


% VÄGGTEMPTERTUREN isolerad vägg
L_i=0.1; %m, isoleringens tjocklek
lambda_i=0.037; %W/mK, mineralulls värmeledningsförmåga
U_i=lambda_i/L_i;
U2=1/(1/U1+1/U_i);

a=-sigma;
b=(-kp+U2);
c=(0.2.*T_u.^4+0.8.*T_a.^4).*sigma+abs*I+kp.*T_u-T_i*U2;

T_w=(0:0.01:30)+273; % i C
y_w=a.*T_w.^4+b.*T_w+c; %y=dT/dt, y=0

% sök STEADYSTATE (nollställe t. y_w)
y_w=y_w.^2;
miny=min(y_w);
index_w=find(y_w<=miny);
zero_Tw2=T_w(index_w);



% Flödet ut ur en oisolerad södervägg per kvadr.meter
F1=U1*(T_i - zero_Tw);

% Flödet ut ur en isolerad typ norrvägg per kvadr.meter
F2=U2*(T_i - zero_Tw2);

% Flödet ut ur en genomsnittlig kvadratmeter vägg, med oisolerad södervägg
U3=0.6086;
F3=U3*(T_i - zero_Tw);

% Flödet ut ur en genomsnittlig kvadratmeter vägg, med isolerad södervägg
U4=0.4680;
F4=U4*(T_i - zero_Tw2);

% Flödet ut ur en genomsnittlig kvadratmeter vägg, med isolerad söder- och västervägg
U5=0.4112;
F5=U5*(T_i - zero_Tw2);


% Så lång tid tar det för väggens yta att komma upp i temperatur.
% t_utan=int_solution(T_u, zero_Tw);
% t_mediso=int_solution(T_u, zero_Tw2);

output=[zero_Tw];




