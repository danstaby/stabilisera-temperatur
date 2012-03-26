function vec = LHST(uu, wu, Tu, freeu, freew, freeT)
global TensX TensZ divxx divzz divzx divxz A QT LM alpha beta nu rho penalty Gu Gw GT


if(nargin == 3)
  freeu = 1:size(uu,1);
  freew = 1:size(wu,1);
  freeT = 1:size(Tu,1);
end

vec = double(ttv( ttv( TensX(freeT,freeu,freeT),Tu,3 ),uu,2) ...
    + ttv( ttv( TensZ(freeT,freew,freeT),Tu,3 ),wu,2)) ...
    + (alpha*A(freeT,freeT)+QT(freeT,freeT))*Tu;