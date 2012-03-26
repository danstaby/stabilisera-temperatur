function theta = angletheta(alpha, beta, gamma)

% gamma = vinkel mellan norr och (fönstrets normal - 90)
beta = beta - 90 - gamma;
if beta<-90 || beta>90 || alpha <0
    theta = 90;
else
    theta = acosd(cosd(beta)*cosd(alpha));
end