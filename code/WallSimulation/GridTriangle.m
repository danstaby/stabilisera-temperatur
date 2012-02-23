classdef GridTriangle < handle
   
    properties
        
        vNodes
        hNodes
        bDraw
    end
    
    
    methods
        function obj = GridTriangle()
            
           obj.bDraw = 1;
        end
        
        function  setNodes(obj, vIDs, vHandles)
            obj.vNodes = vIDs;
            obj.hNodes = vHandles;
        end
        
        function bool = containsNode(Identifier)
            bool = 0;
           for n = 1:3
              if(obj.vNodes(n) == Identifier)
                  bool = 1;
                  break
              end
           end
        end
        
        
        function bool = isInside(obj, pos)
            
            if(obj.bDraw == 0)
                bool = 0;
                return;
            end
            
            X = zeros(3,3);
            p = [pos 0];
            for n=1:3
                
               X(n,:) = [obj.hNodes(n).vPosition, 0];
            end
            
             a = X(1,:);
             b = X(2,:);
             c = X(3,:);
             
             if(GridTriangle.isOnSameSide(p,a,b,c) == 1 && ...
                  GridTriangle.isOnSameSide(p,b,a,c) == 1 && ...
                  GridTriangle.isOnSameSide(p,c,a,b) == 1)
                    bool = 1;
             else
                 bool = 0;
             end
        end
    end
    
    methods(Static)
       function bool = isOnSameSide(p1, p2, a, b)
            cp1 = cross(b-a, p1-a);
            cp2 = cross(b-a, p2-a);
            
            if(dot(cp1, cp2) >= 0)
                bool = 1;
            else
                bool = 0;
            end
        end 
        
    end
    
end