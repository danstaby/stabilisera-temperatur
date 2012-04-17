function U = calcuvalue(L, K)
%U = calcuvalue(L,K)
%
% This function calculates the U Value for the wall specified by
% the column vector L (length) and the column vector K (thermal
% conductivity).


R = sum(L./K);
U = 1/R;

