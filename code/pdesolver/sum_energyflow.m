clc
clear all

% Producerar grafer över det totala energiflödet genom byggnaden över ett
% dygn.
% S+V-väggar, norrväggen, burspråket, 
% taket, solen, grunden, konstanta energiflöden => 8 st
% UT positivt - IN negativt

load rawwalldata.mat
load roofsimulations.mat
    % nosunapril
    % nosundec
    % sunaprilnorth
    % sunaprilsouth
    % sundecnorth
    % sundecsouth
    
% TIDEN
length=length(aprnosun1);
tiden=aprnosun2(length-172:length,1);
tiden=tiden/3600-24*29;

% SOL in gm fönster
load sunpower.mat
    aprsun5mod=-sunpowerapril(2:174,2);
    decsun5mod=-sunpowerapril(2:174,2);

% GRUNDEN och KONSTANT
    foundationdec=3461;
    foundationapr=3575;
    constant=-10500;
    ettor=ones(173,1);
    apr6mod=ettor*(foundationapr+constant);
    dec6mod=ettor*(foundationdec+constant);

A1=151+61; % Area av söder + väster vägg
A2=290; % Area av norrväggen
A3=47; % Area burspråk
A4=257; % Area taket
A5=109; % area av fönster genom vilka solen skiner


% APRIL MOLN
% aprnosun1 oisolerad (s+v)
aprnosun1mod=aprnosun1(length-172:length,2);

% aprnosun2 isolerad (norr)
aprnosun2mod=aprnosun2(length-172:length,2);

% aprnosun3 burspråket
aprnosun3mod=aprnosun3(length-172:length,2);

% nosunapril taket
aprnosun4mod=nosunapril(2:174,2);

% INGEN SOL

% grunden och konstant
apr6mod;

aprnosun_mod=A1*aprnosun1mod+A2*aprnosun2mod+...
    A3*aprnosun3mod+A4*aprnosun4mod+...
    apr6mod;
aprnosun_modi=(A2+A1)*aprnosun2mod+...
    A3*aprnosun3mod+A4*aprnosun4mod+...
    apr6mod;

figure(1)
plot(tiden, aprnosun_mod)
hold on
plot(tiden, aprnosun_modi, '-k')
title('April utan sol. summa energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
hold off


% APRIL SOL
% aprsun1 oisolerad (s+v)
aprsun1mod=aprsun1(length-172:length,2);

% aprsun2 isolerad (norr)
aprsun2mod=aprsun2(length-172:length,2);

% aprsun3 burspråket
aprsun3mod=aprsun3(length-172:length,2);

% nosunapril taket
aprsun4smod=sunaprilsouth(2:174,2);
aprsun4nmod=sunaprilnorth(2:174,2);

% SOL in gm fönster
aprsun5mod;

% GRUND och KONSTANT
apr6mod;

aprsun_mod=A1*aprsun1mod+A2*aprsun2mod+...
    A3*aprsun3mod+A4/2*aprsun4nmod+A4/2*aprsun4smod...
    +A5*aprsun5mod+apr6mod;
aprsun_modi(:,2)=(A1+A2)*aprsun2mod+...
    A3*aprsun3mod+A4/2*aprsun4nmod+A4/2*aprsun4smod+...
    A5*aprsun5mod+apr6mod;

figure(2)
plot(tiden, aprsun_mod, '-r')
hold on
plot(tiden, aprsun_modi, '-k')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
title('April med sol. summa energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)
hold off


% DECEMBER MOLN
% decnosun1 oisolerad (s+v);
decnosun1mod=decnosun1(length-172:length,2);

% decnosun2 isolerad (norr)
decnosun2mod=decnosun2(length-172:length,2);

% decnosun3 burspråket
decnosun3mod=decnosun3(length-172:length,2);

% nosunapril taket
decnosun4mod=nosundec(2:174,2);

% INGEN SOL

% GRUND och KONSTANT
dec6mod;

decnosun_mod=A1*decnosun1mod+A2*decnosun2mod+...
    A3*decnosun3mod+A4*decnosun4mod+...
    dec6mod;
decnosun_modi=(A1+A2)*decnosun2mod+...
    A3*decnosun3mod+A4*decnosun4mod+...
    dec6mod;
figure(3)
plot(tiden, decnosun_mod, '-g')
hold on
plot(tiden, decnosun_modi, '-k')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
title('December utan sol. summa energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)
hold off

% DECEMBER SOL fig 4
% decsun1 oisolerad (s+v)
decsun1mod=decsun1(length-172:length,2);

% decsun2 isolerad (norr)
decsun2mod=decsun2(length-172:length,2);

% decnosun3 burspråket
decsun3mod=decsun3(length-172:length,2);

% decnosun4 taket
decsun4smod=sundecsouth(2:174,2);
decsun4nmod=sundecnorth(2:174,2);

% SOL in gm fönster
decsun5mod;

% GRUND och KONSTANT
dec6mod;

decsun_mod=A1*decsun1mod+A2*decsun2mod+...
    A3*decsun3mod+A4/2*decsun4nmod+A4/2*decsun4smod+...
    A5*decsun5mod+dec6mod;
decsun_modi=(A1+A2)*decsun2mod+...
    A3*decsun3mod+A4/2*decsun4nmod+A4/2*decsun4smod+...
    A5*decsun5mod+dec6mod;
figure(4)
plot(tiden, decsun_mod, '-m')
hold on
plot(tiden, decsun_modi, '-k')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
title('December med sol. summa energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)
hold off