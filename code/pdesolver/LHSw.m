function vec = LHSw(uu, wu, Tu, freeu, freew, freeT)
global TensX TensZ divxx divzz divzx divxz A QT LM alpha beta nu rho penalty Gu Gw GT


if(nargin == 3)
  freeu = 1:size(uu,1);
  freew = 1:size(wu,1);
  freeT = 1:size(Tu,1);
end

vec = double(ttv( ttv( TensX(freew,freeu,freew),wu,3 ),uu,2)) ...
    + double(ttv( ttv( TensZ(freew,freew,freew),wu,3 ),wu,2)) ...
    + nu*A(freew,freew)*wu + penalty*( divzx(freew,freeu)*uu ...
		                 +divzz(freew,freew)*wu )/rho ...
    + LM(freew, freeT)*Tu;