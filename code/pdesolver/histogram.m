clc
clear all

% S+V-väggar, S+V-väggar isolerade, norrväggen, burspråket, 
% taket, solen, grunden, konstanta energiflöden => 8 st

load rawwalldata.mat
load roofsimulations.mat
    % nosunapril
    % nosundec
    % sunaprilnorth
    % sunaprilsouth
    % sundecnorth
    % sundecsouth
load sunpower.mat

colors=colormap(summer(8));

sunpowerapr5mod=-sunpowerapril(2:174);
sunpowerdec5mod=-sunpowerapril(2:174);


% APRIL MOLN
% aprnosun1 oisolerad (s+v)
length=length(aprnosun1);
aprnosun1mod(:,1)=aprnosun1(length-173:length,1);
aprnosun1mod(:,1)=aprnosun1mod(:,1)/3600-24*29;
aprnosun1mod(:,2)=aprnosun1(length-173:length,2);
A1=151+61;

% aprnosun2 isolerad (norr)
aprnosun2mod(:,1)=aprnosun2(length-173:length,1);
aprnosun2mod(:,1)=aprnosun2mod(:,1)/3600-24*29;
aprnosun2mod(:,2)=aprnosun2(length-173:length,2);
A2=290;

% aprnosun3 burspråket
aprnosun3mod(:,1)=aprnosun3(length-173:length,1);
aprnosun3mod(:,1)=aprnosun3mod(:,1)/3600-24*29;
aprnosun3mod(:,2)=aprnosun3(length-173:length,2);
A3=47;

% nosunapril taket
aprnosun4mod=nosunapril(2:173);
A4=257;


aprnosun_mod1(:,2)=A1*aprnosun1mod(:,2)+A2*aprnosun2mod(:,2)+...
    A3*aprnosun3mod(:,2)+A4*aprnosun4mod(:,2);
aprnosun_mod1i(:,2)=(A1+A2)*aprnosun2mod(:,2)+...
    A3*aprnosun3mod(:,2)+A4*aprnosun4mod(:,2);
aprnosun_mod2(:,2)=A2*aprnosun2mod(:,2)+A3*aprnosun3mod(:,2)+...
    A4*aprnosun4mod(:,2);
aprnosun_mod3(:,2)=A2*aprnosun2mod(:,2)+A3*aprnosun3mod(:,2);
aprnosun_mod4(:,2)=A2*aprnosun2mod(:,2);


figure(1)
area(aprnosun1mod(:,1), aprnosun_mod1(:,2),'FaceColor',colors(8,:))
hold on % söder och västerväggarna
area(aprnosun1mod(:,1), aprnosun_mod1i(:,2),'FaceColor',colors(7,:))
hold on % söder och västerväggarna, isolerade
area(aprnosun1mod(:,1), aprnosun_mod2(:,2),'FaceColor',colors(6,:))
hold on % norrväggen
area(aprnosun1mod(:,1), aprnosun_mod3(:,2),'FaceColor',colors(5,:))
hold on % burspråket
area(aprnosun1mod(:,1), aprnosun_mod4(:,2),'FaceColor',colors(4,:))
            % taket
title('April utan sol. summa vägg och tak','FontSize',14)
Legend('Söder- och västväggarna','Söder- och västväggarna isolerade','Tak', 'Burspråk','Norrväggen', 'location', 'SW')
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde W','FontSize',12)
% axis([0 24 0 3000])
hold on


	% DECEMBER MOLN
% decnosun1 oisolerad (s+v)
decnosun1mod(:,1)=decnosun1(length-173:length,1);
decnosun1mod(:,1)=decnosun1mod(:,1)/3600-24*29;
decnosun1mod(:,2)=decnosun1(length-173:length,2);
A1=151+61;

% decnosun2 isolerad (norr)
decnosun2mod(:,1)=decnosun2(length-173:length,1);
decnosun2mod(:,1)=decnosun2mod(:,1)/3600-24*29;
decnosun2mod(:,2)=decnosun2(length-173:length,2);
A2=290;

% decnosun3 burspråket
decnosun3mod(:,1)=decnosun3(length-173:length,1);
decnosun3mod(:,1)=decnosun3mod(:,1)/3600-24*29;
decnosun3mod(:,2)=decnosun3(length-173:length,2);
A3=47;

% nosundecil taket
decnosun4mod=nosundec(2:173);
A4=257;

decnosun_mod1(:,2)=A1*decnosun1mod(:,2)+A2*decnosun2mod(:,2)...
    +A3*decnosun3mod(:,2)+A4*decnosun4mod(:,2);
decnosun_mod1i(:,2)=(A1+A2)*decnosun2mod(:,2)+...
    A3*decnosun3mod(:,2)+A4*decnosun4mod(:,2);
decnosun_mod2(:,2)=A2*decnosun2mod(:,2)+A3*decnosun3mod(:,2)+...
    A4*decnosun4mod(:,2);
decnosun_mod3(:,2)=A2*decnosun2mod(:,2)+A3*decnosun3mod(:,2);
decnosun_mod4(:,2)=A2*decnosun2mod(:,2);

figure(2)
area(decnosun1mod(:,1), decnosun_mod1(:,2),'FaceColor',[.5 .9 .8])
hold on
area(decnosun1mod(:,1), decnosun_mod1i(:,2),'FaceColor',[.5 .9 .6])
hold on
area(decnosun1mod(:,1), decnosun_mod2(:,2),'FaceColor','r')
hold on
area(decnosun1mod(:,1), decnosun_mod3(:,2),'FaceColor','m')
hold on
area(decnosun1mod(:,1), decnosun_mod4(:,2))
title('December utan sol. summa vägg och tak','FontSize',14)
Legend('Söder- och västväggarna','Söder- och västväggarna isolerade','Tak', 'Burspråk','Norrväggen', 'location', 'SW')
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde, W','FontSize',12)
% axis([0 24 0 11000])
hold on


% DECEMBER SOL
% decsun1 oisolerad (s+v)
decsun1mod(:,1)=decsun1(length-172:length,1);
decsun1mod(:,1)=decsun1mod(:,1)/3600-24*29;
decsun1mod(:,2)=decsun1(length-172:length,2);
A1=151+61;

% decsun2 isolerad (norr)
decsun2mod(:,1)=decsun2(length-172:length,1);
decsun2mod(:,1)=decsun2mod(:,1)/3600-24*29;
decsun2mod(:,2)=decsun2(length-172:length,2);
A2=290;

% decnosun3 burspråket
decsun3mod(:,1)=decsun3(length-172:length,1);
decsun3mod(:,1)=decsun3mod(:,1)/3600-24*29;
decsun3mod(:,2)=decsun3(length-172:length,2);
A3=47;

% taket
decsun4smod=sundecsouth(2:172);
decsun4nmod=sundecnorth(2:172);
A4=257;

decsun_mod1(:,2)=A1*decsun1mod(:,2)+A2*decsun2mod(:,2)+...
    A3*decsun3mod(:,2)+(A4/2)*decsun4nmod(:,2)+(A4/2)*decsun4smod(:,2);
decsun_mod1i(:,2)=(A1+A2)*decsun2mod(:,2)+A3*decsun3mod(:,2)+...
    (A4/2)*decsun4nmod(:,2)+(A4/2)*decsun4smod(:,2);
decsun_mod2(:,2)=A2*decsun2mod(:,2)+A3*decsun3mod(:,2)+...
    (A4/2)*decsun4nmod(:,2)+(A4/2)*decsun4smod(:,2);
decsun_mod3(:,2)=A2*decsun2mod(:,2)+A3*decsun3mod(:,2);
decsun_mod4(:,2)=A2*decsun2mod(:,2);

figure(3)
area(decsun1mod(:,1), decsun_mod1(:,2),'FaceColor',[.6 .9 .9])
hold on
area(decsun1mod(:,1), decsun_mod1i(:,2),'FaceColor',[.5 .9 .6])
hold on
area(decsun1mod(:,1), decsun_mod2(:,2),'FaceColor','r')
hold on
area(decsun1mod(:,1), decsun_mod3(:,2),'FaceColor','m')
hold on
area(decsun1mod(:,1), decsun_mod4(:,2))
title('December med sol. summa vägg och tak','FontSize',14)
Legend('Söder- och västväggarna','Söder- och västväggarna isolerade','Tak', 'Burspråk','Norrväggen', 'location', 'SW')
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde, W','FontSize',12)
%axis([0 24 0 11000])
hold on

% APRIL SOL
% aprsun1 oisolerad (s+v)
aprsun1mod(:,1)=aprsun1(length-172:length,1);
aprsun1mod(:,1)=aprsun1mod(:,1)/3600-24*29;
aprsun1mod(:,2)=aprsun1(length-172:length,2);
A1=151+61;

% aprsun2 isolerad (norr)
aprsun2mod(:,1)=aprsun2(length-172:length,1);
aprsun2mod(:,1)=aprsun2mod(:,1)/3600-24*29;
aprsun2mod(:,2)=aprsun2(length-172:length,2);
A2=290;

% aprnosun3 burspråket
aprsun3mod(:,1)=aprsun3(length-172:length,1);
aprsun3mod(:,1)=aprsun3mod(:,1)/3600-24*29;
aprsun3mod(:,2)=aprsun3(length-172:length,2);
A3=47;

% taket
aprsun4smod=sunaprilsouth(2:172);
aprsun4nmod=sunaprilnorth(2:172);
A4=257;

aprsun_mod1(:,2)=A1*aprsun1mod(:,2)+A2*aprsun2mod(:,2)+...
    A3*aprsun3mod(:,2)+(A4/2)*aprsun4nmod(:,2)+(A4/2)*aprsun4smod(:,2);
aprsun_mod1i(:,2)=(A1+A2)*aprsun2mod(:,2)+A3*aprsun3mod(:,2)+...
    (A4/2)*aprsun4nmod(:,2)+(A4/2)*aprsun4smod(:,2);
aprsun_mod2(:,2)=A2*aprsun2mod(:,2)+A3*aprsun3mod(:,2)+...
    (A4/2)*aprsun4nmod(:,2)+(A4/2)*aprsun4smod(:,2);
aprsun_mod3(:,2)=A2*aprsun2mod(:,2)+A3*aprsun3mod(:,2);
aprsun_mod4(:,2)=A2*aprsun2mod(:,2);

figure(4)
area(aprsun1mod(:,1), aprsun_mod1(:,2),'FaceColor',[.6 .9 .9])
hold on
area(aprsun1mod(:,1), aprsun_mod1i(:,2),'FaceColor',[.5 .9 .6])
hold on
area(aprsun1mod(:,1), aprsun_mod2(:,2),'FaceColor','r')
hold on
area(aprsun1mod(:,1), aprsun_mod3(:,2),'FaceColor','m')
hold on
area(aprsun1mod(:,1), aprsun_mod4(:,2))
title('April med sol. summa vägg och tak','FontSize',14)
Legend('Söder- och västväggarna','Söder- och västväggarna isolerade','Tak', 'Burspråk','Norrväggen', 'location', 'SW')
xlabel('Tid, h','FontSize',12)
ylabel('Energiutflöde, W','FontSize',12)
%axis([0 24 0 11000])
hold off

