function ret = triderivative(pts)
%Calculates the derivatives of linear triangular elements.
%For notes of the algorithm see Alberty et. al (1999)

G = [1,1,1;pts];

ret = inv(G)*[0,0;1,0;0,1];