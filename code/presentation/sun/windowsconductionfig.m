function [time,qapr,qdec]=windowsconductionfig()
Time = [6;16;30;40;54;64];
Temp = [6;9;6;9;6;9];
t = 6:0.1:64;
T = interp1(Time,Temp,t);
T = T(18*10:42*10);
t = 0:0.1:24;

U = 1;
h = 6.19;
Ueff = 1/(1/(U) + 1/(h));
Tin = 20;

q = Ueff.*(20 - T);
qapr=q;

hold off
fig4 = figure(4);
plot(t, q)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
title('Kyleffekt, värmeledning fönster', 'fontsize', 12)
xlim([0 24])
set(fig4, 'Position', [300 700 300 200])

Time = [6;16;30;40;54;64];
Temp = [-11;-5;-11;-5;-11;-5];
t = 6:0.1:64;
T = interp1(Time,Temp,t);
T = T(18*10:42*10);
t = 0:0.1:24;
h = 35;
Ueff = 1/(1/(U) + 1/(h));
q = Ueff.*(20 - T);

fig5 = figure(5);
plot(t,q)
title('Kyleffekt,  värmeledning fönster', 'fontsize', 12)
xlabel('Tid i timmar')
ylabel('Kyleffekt (W m^{-2})')
xlim([0 24])
set(fig5, 'Position', [800 700 300 200])

qdec=q;
time=t;