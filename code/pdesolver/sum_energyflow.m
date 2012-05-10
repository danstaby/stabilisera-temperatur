clc
clear all

% Producerar grafer över det totala energiflödet genom byggnaden över ett
% dygn.
% S+V-väggar, norrväggen, burspråket, 
% taket, solen, grunden, konstanta energiflöden => 8 st
% UT positivt - IN negativt

colors=colormap(hsv(9));
colors2=colormap(hsv(9));

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
    aprsun=-sunpowerapril(2:174,2);
    decsun5mod=-sunpowerapril(2:174,2);
    decsun=-sunpowerapril(2:174,2);

% strålning UT gm fönaster.
sigma=5.67*10^(-8);
intemp=ones(173,1)*(273+20);
utetemp=temperature();
tempapr=utetemp(:,1)+273;
tempdec=utetemp(:,2)+273;
utapr=0.9*sigma.*(intemp.^4-tempapr.^4);
utdec=0.9*sigma.*(intemp.^4-tempdec.^4);

% GRUNDEN och KONSTANT
    foundationdec=3461;
    foundationapr=3575;
    constant=-10500;
    ettor=ones(173,1);
    apr6mod=ettor*(foundationapr+constant);
    dec6mod=ettor*(foundationdec+constant);
    aprFoundation=ettor*foundationapr;
    decFoundation=ettor*foundationdec;
    constant=ettor*constant;
    

A1=151+61; % Area av söder + väster vägg
A2=290; % Area av norrväggen
A3=47; % Area burspråk
A4=257; % Area taket
A5=109; % area av fönster genom vilka solen skiner
A6=109+89+9; % area alla fönster.



% APRIL MOLN
% aprnosun1 oisolerad (s+v)
swWall=A1*aprnosun1(length-172:length,2);

% aprnosun2 isolerad (norr)
nWall=A2*aprnosun2(length-172:length,2);
swWalli=A1*aprnosun2(length-172:length,2);

% aprnosun3 burspråket
bWall=A3*aprnosun3(length-172:length,2);

% nosunapril taket
roof=A4*nosunapril(2:174,2);

% fönster, INGEN SOL
window=utdec*A6;

% grunden och konstant
foundation=aprFoundation;
constant;

% Positiva energiflöden (UT)
% norrväggen, söder- och västerväggen, burspråket, taket, grunden, fönstern
summaUt1=foundation+window+roof+bWall+nWall+swWall;
summaUt2=foundation+window+roof+bWall+nWall+swWalli;
summaUt3=foundation+window+roof+bWall+nWall;
summaUt4=foundation+window+roof+bWall;
summaUt5=foundation+window+roof;
summaUt6=foundation+window;
summaUt7=foundation;

% Negativa energiflöden (IN)
% solinstrålning, konstanta värmekällor
summaIn1=constant;

% Totalt
total=summaUt1+summaIn1;

figure(1)
area(tiden, summaUt1,'FaceColor',colors(1,:));
hold on
area(tiden, summaUt2,'FaceColor',colors(2,:));
area(tiden, summaUt3,'FaceColor',colors(3,:));
area(tiden, summaUt4,'FaceColor',colors(4,:));
area(tiden, summaUt5,'FaceColor',colors(5,:));
area(tiden, summaUt6,'FaceColor',colors(6,:));
area(tiden, summaUt7,'FaceColor',colors(7,:));
area(tiden, summaIn1,'FaceColor',colors(8,:));
plot(tiden, total, '-k','linewidth', 3);

Legend('Som idag','Med isolering',...
    'Norrväggen', 'Burspårket', 'Taket', 'Fönster',...
    'Grunden', 'Konstant','location', 'SE')
title('April utan sol. Summering av energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)

hold off





% APRIL SOL
% aprsun1 oisolerad (s+v)
swWall=A1*aprsun1(length-172:length,2);

swWallPos=swWall;
swWallNeg=swWall;
for i=1:173
    if swWall(i)<0
        swWallPos(i)=0;
    else
        swWallNeg(i)=0;
    end
end


% aprsun2 isolerad (norr)
nWall=A2*aprsun2(length-172:length,2);
swWalli=A1*aprsun2(length-172:length,2);

swWalliPos=swWalli;
swWalliNeg=swWalli;
for i=1:173
    if swWalli(i)<0
        swWalliPos(i)=0;
    else
        swWalliNeg(i)=0;
    end
end

nWallPos=nWall;
nWallNeg=nWall;
for i=1:173
    if nWall(i)<0
        nWallPos(i)=0;
    else
        nWallNeg(i)=0;
    end
end


% aprsun3 burspråket
bWall=A3*aprsun3(length-172:length,2);
bWallPos=bWall;
bWallNeg=bWall;
for i=1:173
    if bWall(i)<0
        bWallPos(i)=0;
    else
        bWallNeg(i)=0;
    end
end


% nosunapril taket
roof=A4/2*sunaprilsouth(2:174,2)+A4/2*sunaprilnorth(2:174,2);
roofPos=roof;
roofNeg=roof;
for i=1:173
    if roof(i)<0
        roofPos(i)=0;
    else
        roofNeg(i)=0;
    end
end

% fönster, sol in och strålning ut
sunIn=A5*aprsun;
windowOut=A6*utapr;
window=windowOut+sunIn;
windowPos=window;
windowNeg=window;
for i=1:173
    if window(i)<0
        windowPos(i)=0;
    else
        windowNeg(i)=0;
    end
end

% GRUND och KONSTANT
foundation=aprFoundation;
constant;

% Positiva energiflöden (UT)
% norrväggen, söder- och västerväggen, burspråket, taket, grunden, fönstern
summaUt1=foundation+windowPos+roofPos+bWallPos+nWallPos+swWallPos;
summaUt2=foundation+windowPos+roofPos+bWallPos+nWallPos+swWalliPos;
summaUt3=foundation+windowPos+roofPos+bWallPos+nWallPos;
summaUt4=foundation+windowPos+roofPos+bWallPos;
summaUt5=foundation+windowPos+roofPos;
summaUt6=foundation+windowPos;
summaUt7=foundation;

% Negativa energiflöden (IN)
% solinstrålning, konstanta värmekällor
summaIn1=constant+windowNeg+roofNeg+bWallNeg+nWallNeg+swWallNeg;
summaIn2=constant+windowNeg+roofNeg+bWallNeg+nWallNeg+swWalliNeg;
summaIn3=constant+windowNeg+roofNeg+bWallNeg+nWallNeg;
summaIn4=constant+windowNeg+roofNeg+bWallNeg;
summaIn5=constant+windowNeg+roofNeg;
summaIn6=constant+windowNeg;
summaIn7=windowNeg;

% Totalt
total=summaUt1+summaIn1;

figure(2)
area(tiden, summaUt1,'FaceColor',colors(1,:));
hold on
area(tiden, summaUt2,'FaceColor',colors(2,:));
area(tiden, summaUt3,'FaceColor',colors(3,:));
area(tiden, summaUt4,'FaceColor',colors(4,:));
area(tiden, summaUt5,'FaceColor',colors(5,:));
area(tiden, summaUt6,'FaceColor',colors(6,:));
area(tiden, summaUt7,'FaceColor',colors(7,:)); % Grunden

area(tiden, summaIn7,'FaceColor',colors2(8,:)); % Konstanter
area(tiden, summaIn1,'FaceColor',colors2(1,:));
area(tiden, summaIn2,'FaceColor',colors2(2,:));
area(tiden, summaIn3,'FaceColor',colors2(3,:));
area(tiden, summaIn4,'FaceColor',colors2(4,:));
area(tiden, summaIn5,'FaceColor',colors2(5,:));
area(tiden, summaIn6,'FaceColor',colors2(8,:));
area(tiden, summaIn7,'FaceColor',colors2(6,:)); % Konstanter
plot(tiden, total, '-k','linewidth', 3);

Legend('Som idag','Med isolering',...
    'Norrväggen', 'Burspårket', 'Taket', 'Fönster',...
    'Grunden', 'Konstant','location', 'SE')
title('April med sol. Summering av energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)

hold off





% DECEMBER MOLN
% decnosun1 oisolerad (s+v);
swWall=A1*decnosun1(length-172:length,2);

% decnosun2 isolerad (norr)
decnosun2mod=decnosun2(length-172:length,2);
nWall=A2*decnosun2mod;
swWalli=A1*decnosun2mod;

% decnosun3 burspråket
decnosun3mod=decnosun3(length-172:length,2);
bWall=A3*decnosun3mod;

% nosunapril taket
decnosun4mod=nosundec(2:174,2);
roof=A4*decnosun4mod;

% fönster, INGEN SOL
window=utdec*A6;

% GRUND och KONSTANT
foundation=decFoundation;
const=dec6mod;

% Positiva energiflöden (UT)
% norrväggen, söder- och västerväggen, burspråket, taket, grunden, fönstern
summaUt1=foundation+window+roof+bWall+nWall+swWall;
summaUt2=foundation+window+roof+bWall+nWall+swWalli;
summaUt3=foundation+window+roof+bWall+nWall;
summaUt4=foundation+window+roof+bWall;
summaUt5=foundation+window+roof;
summaUt6=foundation+window;
summaUt7=foundation;

% Negativa energiflöden (IN)
% solinstrålning, konstanta värmekällor
summaIn1=constant;

% Totalt
total=summaUt1+summaIn1;

figure(3)
area(tiden, summaUt1,'FaceColor',colors(1,:));
hold on
area(tiden, summaUt2,'FaceColor',colors(2,:));
area(tiden, summaUt3,'FaceColor',colors(3,:));
area(tiden, summaUt4,'FaceColor',colors(4,:));
area(tiden, summaUt5,'FaceColor',colors(5,:));
area(tiden, summaUt6,'FaceColor',colors(6,:));
area(tiden, summaUt7,'FaceColor',colors(7,:));
area(tiden, summaIn1,'FaceColor',colors(8,:));
plot(tiden, total, '-k','linewidth', 3);

Legend('Som idag','Med isolering',...
    'Norrväggen', 'Burspårket', 'Taket', 'Fönster',...
    'Grunden', 'Konstant','location', 'SE')
title('December utan sol. Summering av energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)

hold off




% DECEMBER SOL fig 4
% decsun1 oisolerad (s+v)
decsun1mod=decsun1(length-172:length,2);
swWall=A1*decsun1mod;

% decsun2 isolerad (norr)
decsun2mod=decsun2(length-172:length,2);
nWall=A2*decsun2mod;
swWalli=A1*decsun2mod;

% decnosun3 burspråket
decsun3mod=decsun3(length-172:length,2);
bWall=A3*decsun3mod;

% decnosun4 taket
decsun4smod=sundecsouth(2:174,2);
decsun4nmod=sundecnorth(2:174,2);
roof=A4*decsun4nmod/2+A4*decsun4smod/2;

% fönster, sol in och strålning ut
sunIn=A5*decsun;
windowOut=A6*utdec;
window=windowOut+sunIn;
windowPos=window;
windowNeg=window;
for i=1:173
    if window(i)<0
        windowPos(i)=0;
    else
        windowNeg(i)=0;
    end
end

% GRUND och KONSTANT
foundation=decFoundation;
constant;

% Positiva energiflöden (UT)
% norrväggen, söder- och västerväggen, burspråket, taket, grunden, fönstern
summaUt1=foundation+windowPos+roof+bWall+nWall+swWall;
summaUt2=foundation+windowPos+roof+bWall+nWall+swWalli;
summaUt3=foundation+windowPos+roof+bWall+nWall;
summaUt4=foundation+windowPos+roof+bWall;
summaUt5=foundation+windowPos+roof;
summaUt6=foundation+windowPos;
summaUt7=foundation;

% Negativa energiflöden (IN)
% solinstrålning, konstanta värmekällor
summaIn1=constant+windowNeg;
summaIn2=windowNeg;

% Totalt
total=summaUt1+summaIn1;

figure(4)
area(tiden, summaUt1,'FaceColor',colors(1,:));
hold on
area(tiden, summaUt2,'FaceColor',colors(2,:));
area(tiden, summaUt3,'FaceColor',colors(3,:));
area(tiden, summaUt4,'FaceColor',colors(4,:));
area(tiden, summaUt5,'FaceColor',colors(5,:));
area(tiden, summaUt6,'FaceColor',colors(6,:));
area(tiden, summaUt7,'FaceColor',colors(7,:));
area(tiden, summaIn1,'FaceColor',colors(8,:));
area(tiden, summaIn2,'FaceColor',colors(6,:));
plot(tiden, total, '-k','linewidth', 3);

Legend('Som idag','Med isolering',...
    'Norrväggen', 'Burspårket', 'Taket', 'Fönster',...
    'Grunden', 'Konstant','location', 'NE')
title('December med sol. Summering av energiflöden','FontSize',14)
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)

hold off
