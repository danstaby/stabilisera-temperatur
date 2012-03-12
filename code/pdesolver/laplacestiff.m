function M = laplacestiff(vertices)
%Calculates the local stiffness matrix. This Function
%is described by Alberty et. al. (1999) "Remarks around 50 lines of
%Matlab: short finite element implementation"


G = [ones(1,3);vertices] \ [zeros(1,2);eye(2)];
M = det([ones(1,3);vertices]) * G * G'/2;