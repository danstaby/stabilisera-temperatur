function ret = zakianmethod()
global x

a=[(1.283767675e+1)+i*1.666063445,(1.222613209e+1)+i*5.012718792,...
    (1.09343031e+1)+i*8.40967312,8.77643472+i*1.19218539e+1,...
    5.22545336+i*1.57295290e+1];

K=[(-3.69020821e+4)+i*1.96990426e+5,(6.12770252e+4)-i*9.54086255e+4,...
      -(2.89165629e+4)+i*1.81691853e+4,(4.65536114e+3)-i*1.90152864,...
      -(1.18741401e+2)-i*1.41303691e+2];

x = [0:0.01:0.5]';

figure(3)
plot(x, 2*x*20)
hold on

for t = 1000:1000:30000
  
  fval = zeros(size(x,1),1);

  for n = 1:5
    fval = fval + real(K(n)*Fs(a(n)/t));
  end
  
  fval = 2*fval/t;

  plot(x,fval, 'r')
  if(t == 1000)
    plot(x,fval, 'm')
  end
end

plot(x,fval, 'k')
hold off

function val = Fs(s)
global x
Tu = 10;
T0 = 0;
Ti = 20;
L = 0.5;
alpha = 5.2e-7;

val = (Tu-T0)*sinh((L-x)*sqrt(s/alpha))/(s*sinh(L*sqrt(s/alpha))) + ...
      (Ti-T0)*x/(s*L) + T0/s;


function fun = func(t)

fun = t.^2+2;