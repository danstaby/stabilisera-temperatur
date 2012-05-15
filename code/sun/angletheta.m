function theta = angletheta(alpha, beta, gamma)
% gamma: the angle between north and the normal of the window
% öster om norr positivt

beta = beta - gamma;

if beta<-90 || beta>90 || alpha <0
    theta = 90;
else
    theta = acosd(cosd(beta)*cosd(alpha));
end