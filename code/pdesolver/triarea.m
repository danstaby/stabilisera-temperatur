function area = triarea(pts)
%This function calculates the area of a triangle

M = [1,1,1;pts];
area = 0.5*det(M);

