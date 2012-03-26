function vec = LHSu(uu, wu, Tu, freeu, freew, freeT)

global TensX TensZ divxx divzz divzx divxz A QT LM alpha beta nu rho penalty Gu Gw GT
if(nargin == 3)
  freeu = 1:size(uu,1);
  freew = 1:size(wu,1);
  freeT = 1:size(Tu,1);
end


vec = double(ttv( ttv( TensX(freeu,freeu,freeu),uu,3 ),uu,2)) ...
    + double(ttv( ttv( TensZ(freeu,freew,freeu),uu,3 ),wu,2)) ...
    + nu*A(freeu,freeu)*uu + penalty*( divxx(freeu,freeu)*uu ...
    + divxz(freeu,freew)*wu )/rho;
