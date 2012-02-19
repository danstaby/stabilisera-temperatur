classdef LinkedList < handle
    
    properties(SetAccess = private)
        iSize
    end
    
    properties(GetAccess=private)
        hFirstNode
        hLastNode
        iRollingID
        
    end
    
    methods
        function obj = LinkedList()
            obj.iRollingID = 0;
        end
        
        function id = addElement(obj, Value, Identifier)
            hNode = LinkedNode;
            hNode.Value = Value;
            
            if(nargin == 2)
                hNode.iIdentifier = obj.getNewRollingID();
            else
                hNode.iIdentifier = Identifier;
            end
            
            id = hNode.iIdentifier;
            
            if(obj.iSize == 0)
                obj.hFirstNode = hNode;
                obj.hLastNode = hNode;
            else
                obj.hLastNode.hNextNode = hNode;
                hNode.hPrevNode = obj.hLastNode;
                obj.hLastNode = hNode;
            end
            
            obj.iSize = obj.iSize + 1;
        end
        
        function Value = getElement(obj, Identifier)
           hNode = obj.findElement(Identifier);
           
           if(hNode == 0)
               Value = 0;
           else
               Value = hNode.Value;
           end
        end
        
        function Value = getFirstElement(obj)
            if(obj.hFirstNode ~= 0)
                Value = obj.hFirstNode.Value;
            else
                Value = 0;
            end
        end
        
        function [Value ID] = getFirstElementAndID(obj)
            if(obj.hFirstNode ~= 0)
                Value = obj.hFirstNode.Value;
                ID = obj.hFirstNode.iIdentifier;
            else
                Value = 0;
                ID = 0;
            end
        end
        
        function ID = getFirstElementID(obj)
            if(obj.hFirstNode ~= 0)
                ID = obj.hFirstNode.iIdentifier;
            else
                ID = 0;
            end
        end
        
        function ID = getNextElementID(obj, Identifier)
            ID = 0;
            hNode = findElement(Identifier);
            if(hNode ~= 0)
               ID =  hNode.iIdentifier;
            end
        end
        
        function Value = getNextElement(obj, Identifier)
            Value = 0;
            hNode = findElement(Identifier);
            if(hNode ~= 0)
                Value = hNode.Value;
            end
        end
        
        function [Value ID] = getNextElementAndID(obj, Identifier)
            Value = 0;
            ID = 0;
           
            hNode = obj.findElement(Identifier);
            if(hNode ~= 0)
                if(hNode.hNextNode ~= 0)
                    Value = hNode.hNextNode.Value;
                    ID = hNode.hNextNode.iIdentifier;
                end
            end
        end
        
        function clearList(obj)
            obj.iRollingID = 0;
            obj.hFirstNode = 0;
            obj.hLastNode = 0;
            obj.iSize = 0;
        end
        
        function removeElement(obj, Identifier)
            hNode = obj.findElement(Identifier);
            
            if(hNode ~= 0)
                obj.iSize = obj.iSize - 1;
                
                hNext = hNode.hNextNode;
                hPrev = hNode.hPrevNode;
                
                if(hNext == 0)
                    obj.hLastNode = hPrev;
                else
                    hNext.hPrevNode = hPrev;
                end
                
                if(hPrev == 0)
                    obj.hFirstNode = hNext;
                else
                    hPrev.hNextNode = hNext;
                end
            end
        end
    end
    
    methods(Access = private)
        function hNode = findElement(obj, Identifier)
            
            hNode = 0;
            
            hN = obj.hFirstNode;
            
            while(hN ~= 0)
                if(hN.iIdentifier == Identifier)
                   hNode = hN;
                   hN = 0;
                else
                    hN = hN.hNextNode;
                end
            end
        end
        
        function ID = getNewRollingID(obj)
           obj.iRollingID = obj.iRollingID + 1;
           ID = obj.iRollingID;
        end
    end
    
    
end