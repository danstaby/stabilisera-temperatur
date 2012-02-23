classdef GridNode < handle
    
   properties
       vPosition
       bBoundary
       
   end
   
   properties(SetAccess = private)
      cBoundaryConditions
   end
   
   methods
       function obj = GridNode()
          obj.vPosition = 0;
          
          obj.cBoundaryConditions = {};
          
       end
       
       function addBoundaryCondition(obj, strFunc, dOrder, normal)
            s = max(size(obj.cBoundaryConditions));
            obj.cBoundaryConditions{s+1} = {strFunc, dOrder, normal};
       end
       
       function [f fprim] = getBoundaryFunction(obj ,vTangent, iParameter)
           f = 0;
           fprim = 0;
           if(nargin < 3)
               iParameter = 0;
           end
          
          for n = 1:max(size(obj.cBoundaryConditions))
             row = obj.cBoundaryConditions{n};
             normal = row{3};
             
             if(dot(normal,vTangent) == 0)
                 strFun = row{1};
                 strFun = strrep(strFun, 'VALUE', num2str(iParameter));
                 fun = inline(strFun);
                 
                 if(row{2} == 0)
                     f = fun;
                 else
                     fprim = fun;
                 end
             end
              
          end
       end
       
   end
    
end