% % Naturkonstanter för luft


% TERMISK DIFFUSITET
alpha=k/(rho*c_p);

% k är termisk konduktivitetet, W/(mK)
% rho är densiteten, kg/m^3
% c_p specifika värmekapaciteten, J/(kgK)

% Luft: 
alpha_l=1.9*10^(-5); %m/s^2



% TERMISK KONDUKTIVITET (Värmeledningsförmåga), k

% Eftersom beräknat och uppmätt skiljer sig så är det svårt att få en bra
% 'formel'. Den varierar mellan 0.0204 och 0.0299 för -50 till 80 °C

% Vid rumstemperatur: 
% Beräknat: 0,024 W/mK
% uppmätt: 0,031 W/mK, enligt Schroeder 

% Värmeledningsförmågan ändras med temperaturen. 
% För de flesta ämnen minskar den något med stigande temperatur. 
% Vid låga temperaturer kan det även bero på trycket.



% DENSITET
% Ideala gaslagen
% p*V = n*R*T;
molmassa_luft=0.029; % kg
rho=p/(R*T*molmassa_luft);

p; %lufttrycket, N/m2
R = 8.3145; %gaskonstanten,  J/mol K
T; % absoluta temperaturen (i Kelvin).

V; % gasens volym (i m3)
n; % substansmängd eller molantal (i mol)
    


% SPECIFIK VÄRMEKAPACITET
% kan approximeras med konstanten 
c_v=1.005; % kJ/kgK

% Saknar data för hur den beror av fuktigheten.



% VOLYMETRISKA EXPANSIONSKOEFFICIENTEN
% använd allmänna gaslagen (ty fuktighet påverkar så lite) 
beta=dV/dT;
% Allmänna gaslagen, V = n*R*T/p, ger då:
beta=n*R/p;

% Expanderar linjärt (konstant)



% KINEMATISKA VISKOSITETEN

% Den kinematiska viskositeten $\nu$ definieras som: 
nu = mu/rho;

% rho: densiteten
rho=p/(R*T*molmassa_luft);

% Dynamisk viskositet,mu, luft:
% 1,67*10^(-5) Pa s, vid 0°C
% 1,78*10^(-5) Pa s, vid 15°C
% Den dynamiska viskositeten för gasen kan anses vara 
% OBEROENDE AV TRYCK i de flesta tekniska tillämpningar.

% Graf över beroende av tryck och temperatur:
% http://en.wikipedia.org/wiki/File:Air_dry_dynamic_visocity_on_pressure_te
% mperature.svg
% Bedöm gärna själv om det verkar rimligt att bortse ifrån.

% Sutherlands formel kan användas för att bestämma den 
% dynamiska viskositeten för en ideal gas som funktion av temperatur.

% \begin{equation}
mu=mu0*(T_0+C)/(T+C)*(T/T0)^(3/2);

% Luft: 
C=120; %K 
T_0=291.15; %K
mu_0=18.27; % Pa s

% Riktig för temperaturer 0 < T < 555 K 
% med ett fel på grund av tryck på mindre än 10\% under 3.45 MPa.
