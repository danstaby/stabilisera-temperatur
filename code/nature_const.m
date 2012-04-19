% % Naturkonstanter f�r luft


% TERMISK DIFFUSITET
alpha=k/(rho*c_p);

% k �r termisk konduktivitetet, W/(mK)
% rho �r densiteten, kg/m^3
% c_p specifika v�rmekapaciteten, J/(kgK)

% Luft: 
alpha_l=1.9*10^(-5); %m/s^2



% TERMISK KONDUKTIVITET (V�rmeledningsf�rm�ga), k

% Eftersom ber�knat och uppm�tt skiljer sig s� �r det sv�rt att f� en bra
% 'formel'. Den varierar mellan 0.0204 och 0.0299 f�r -50 till 80 �C

% Vid rumstemperatur: 
% Ber�knat: 0,024 W/mK
% uppm�tt: 0,031 W/mK, enligt Schroeder 

% V�rmeledningsf�rm�gan �ndras med temperaturen. 
% F�r de flesta �mnen minskar den n�got med stigande temperatur. 
% Vid l�ga temperaturer kan det �ven bero p� trycket.



% DENSITET
% Ideala gaslagen
% p*V = n*R*T;
molmassa_luft=0.029; % kg
rho=p/(R*T*molmassa_luft);

p; %lufttrycket, N/m2
R = 8.3145; %gaskonstanten,  J/mol K
T; % absoluta temperaturen (i Kelvin).

V; % gasens volym (i m3)
n; % substansm�ngd eller molantal (i mol)
    


% SPECIFIK V�RMEKAPACITET
% kan approximeras med konstanten 
c_v=1.005; % kJ/kgK

% Saknar data f�r hur den beror av fuktigheten.



% VOLYMETRISKA EXPANSIONSKOEFFICIENTEN
% anv�nd allm�nna gaslagen (ty fuktighet p�verkar s� lite) 
beta=dV/dT;
% Allm�nna gaslagen, V = n*R*T/p, ger d�:
beta=n*R/p;

% Expanderar linj�rt (konstant)



% KINEMATISKA VISKOSITETEN

% Den kinematiska viskositeten $\nu$ definieras som: 
nu = mu/rho;

% rho: densiteten
rho=p/(R*T*molmassa_luft);

% Dynamisk viskositet,mu, luft:
% 1,67*10^(-5) Pa s, vid 0�C
% 1,78*10^(-5) Pa s, vid 15�C
% Den dynamiska viskositeten f�r gasen kan anses vara 
% OBEROENDE AV TRYCK i de flesta tekniska till�mpningar.

% Graf �ver beroende av tryck och temperatur:
% http://en.wikipedia.org/wiki/File:Air_dry_dynamic_visocity_on_pressure_te
% mperature.svg
% Bed�m g�rna sj�lv om det verkar rimligt att bortse ifr�n.

% Sutherlands formel kan anv�ndas f�r att best�mma den 
% dynamiska viskositeten f�r en ideal gas som funktion av temperatur.

% \begin{equation}
mu=mu0*(T_0+C)/(T+C)*(T/T0)^(3/2);

% Luft: 
C=120; %K 
T_0=291.15; %K
mu_0=18.27; % Pa s

% Riktig f�r temperaturer 0 < T < 555 K 
% med ett fel p� grund av tryck p� mindre �n 10\% under 3.45 MPa.
