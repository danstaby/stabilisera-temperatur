function theta = angletheta(alpha, beta, gamma)
% gamma = angle of windows normal relative north, east positive

beta = beta - gamma;

if beta<-90 || beta>90 || alpha <0
    theta = 90;
else
    theta = acosd(cosd(beta)*cosd(alpha));
end