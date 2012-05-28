close all
clc

a1 = 1.09;
b1 = 0.23;
c1 = 1;

a2 = 0;
b2 = 0.53;
c2 = 0.78;

v1 = 0:0.01:4.88;
v2 = 4.88:0.01:20;
T = 273.15 + [-20:10:30];
%color = ['b','r','g','m','c','y'];

for i = 1:length(T)
    h1 = 5.678.*(a1 + b1.*(965.42.*v1./T(i)).^c1);
    h2 = 5.678.*(a2 + b2.*(965.42.*v2./T(i)).^c2);
    plot(v1, h1, v2, h2, 'b')
    hold on
end

ymin = 0;
ymax = 90;
xval = 4.88;

plot([xval, xval],[ymin, ymax], 'k:')
xlabel('v, ms^{-1}')
ylabel('h, Wm^{-2}K^{-1}')