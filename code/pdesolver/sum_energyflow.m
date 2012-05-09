clc
clear all

% Producerar grafer över det totala energiflödet genom byggnaden över ett
% dygn.

load rawwalldata.mat
load roofsimulations.mat
    % nosunapril
    % nosundec
    % sunaprilnorth
    % sunaprilsouth
    % sundecnorth
    % sundecsouth

    
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

aprnosun_mod(:,2)=A1*aprnosun1mod(:,2)+A2*aprnosun2mod(:,2)+A3*aprnosun3mod(:,2)+A4*aprnosun4mod(:,2);
aprnosun_modi(:,2)=(A2+A1)*aprnosun2mod(:,2)+A3*aprnosun3mod(:,2)+A4*aprnosun4mod(:,2);

figure(1)
plot(aprnosun1mod(:,1), aprnosun_mod(:,2))
hold on
plot(aprnosun1mod(:,1), aprnosun_modi(:,2), '-k')
title('April utan sol. summa vägg och tak')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
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

% aprsun3 burspråket
aprsun3mod(:,1)=aprsun3(length-172:length,1);
aprsun3mod(:,1)=aprsun3mod(:,1)/3600-24*29;
aprsun3mod(:,2)=aprsun3(length-172:length,2);
A3=47;

% nosunapril taket
aprsun4smod=sunaprilsouth(2:172);
aprsun4nmod=sunaprilnorth(2:172);
A4=257;

aprsun_mod(:,2)=A1*aprsun1mod(:,2)+A2*aprsun2mod(:,2)+A3*aprsun3mod(:,2)+A4/2*aprsun4nmod(:,2)+A4/2*aprsun4smod(:,2);
aprsun_modi(:,2)=(A1+A2)*aprsun2mod(:,2)+A3*aprsun3mod(:,2)+A4/2*aprsun4nmod(:,2)+A4/2*aprsun4smod(:,2);

figure(2)
plot(aprsun1mod(:,1), aprsun_mod(:,2), '-r')
hold on
plot(aprsun1mod(:,1), aprsun_modi(:,2), '-k')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
title('April med sol. summa vägg och tak')



% DECEMBER MOLN
% decnosun1 oisolerad (s+v);
decnosun1mod(:,1)=decnosun1(length-172:length,1);
decnosun1mod(:,1)=decnosun1mod(:,1)/3600-24*29;
decnosun1mod(:,2)=decnosun1(length-172:length,2);
A1=151+61;

% decnosun2 isolerad (norr)
decnosun2mod(:,1)=decnosun2(length-172:length,1);
decnosun2mod(:,1)=decnosun2mod(:,1)/3600-24*29;
decnosun2mod(:,2)=decnosun2(length-172:length,2);
A2=290;

% decnosun3 burspråket
decnosun3mod(:,1)=decnosun3(length-172:length,1);
decnosun3mod(:,1)=decnosun3mod(:,1)/3600-24*29;
decnosun3mod(:,2)=decnosun3(length-172:length,2);
A3=47;

% nosunapril taket
decnosun4mod=nosundec(2:173);
A4=257;

decnosun_mod(:,2)=A1*decnosun1mod(:,2)+A2*decnosun2mod(:,2)+A3*decnosun3mod(:,2)+A4*decnosun4mod(:,2);
decnosun_modi(:,2)=(A1+A2)*decnosun2mod(:,2)+A3*decnosun3mod(:,2)+A4*decnosun4mod(:,2);
figure(3)
plot(decnosun1mod(:,1), decnosun_mod(:,2), '-g')
hold on
plot(decnosun1mod(:,1), decnosun_modi(:,2), '-k')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
title('December utan sol. summa vägg och tak')
hold on



% DECEMBER SOL fig 4
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

decsun4smod=sundecsouth(2:172);
decsun4nmod=sundecnorth(2:172);
A4=257;

decsun_mod(:,2)=A1*decsun1mod(:,2)+A2*decsun2mod(:,2)+A3*decsun3mod(:,2)+A4/2*decsun4nmod(:,2)+A4/2*decsun4smod(:,2);
decsun_modi(:,2)=(A1+A2)*decsun2mod(:,2)+A3*decsun3mod(:,2)+A4/2*decsun4nmod(:,2)+A4/2*decsun4smod(:,2);
figure(4)
plot(decsun1mod(:,1), decsun_mod(:,2), '-m')
hold on
plot(decsun1mod(:,1), decsun_modi(:,2), '-k')
Legend('Som idag', 'Med isolerad norr- och västervägg', 'location', 'SW')
title('December med sol. summa vägg och tak')
