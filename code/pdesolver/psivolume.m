function V = psivolume(P)
%Volume of our test function
%for description of algorithm see http://mathworld.wolfram.com/Tetrahedron.html

M = ones(4,4);
M(1:3, 3) = 0;
M(1:3,1:2) = P';
M(4, 1:3) = [P(1:2,1)' 1];
V = abs(det(M))/6;