function g = gvalue(g0, p, q, theta)
% g0 = g-value at normal incidence
% p and q  depending on window properties
% theta = suns angle relative windows normal

a = 8;
b = 0.25/q;
c = 1 - a - b;

alfa = 5.2 + 0.7*q;
beta = 2;
gamma = 5.26 + 0.06*p + (0.73 + 0.04*p)*q;

z = theta/90;

g = g0*(1 - a*z^(alfa) - b*z^(beta) - c*z^(gamma));
