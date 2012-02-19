classdef GridNode < handle
    
   properties
       vPosition
       bBoundary
   end
   
   properties(SetAccess = private)
       bNeumann
       bDirichlet
       fNeumann
       fDirichlet
   end
   
   methods
       function obj = GridNode()
          obj.bBoundary = 0;
          obj.bNeumann = 0;
          obj.bDirichlet = 0;
          obj.vPosition = 0;
       end
       
       function setNeumann(obj, strFunc)
           obj.bNeumann = 1;
           obj.fNeumann = inline(strFunc);
       end
       
       function setDirichlet(obj, strFunc)
           obj.bDirichlet = 1;
           obj.fDirichlet = inline(strFunc);
       end
       
   end
    
end