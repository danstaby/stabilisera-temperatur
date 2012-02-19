classdef LinkedNode < handle
    
   properties
       hNextNode
       hPrevNode
       Value
       iIdentifier
   end
   
   methods
       function obj = LinkedNode(val)
            if(nargin ~= 0)
                obj.Value = val;
            else
                obj.Value = 0;
            end
            obj.hPrevNode = 0;
            obj.hNextNode = 0;
       end
       
       function val = get.Value(obj)
            val = obj.Value;
       end
       
       function set.Value(obj,val)
            obj.Value = val;
       end
       
       function id = get.iIdentifier(obj)
           id = obj.iIdentifier;
       end
       
       function set.iIdentifier(obj, id)
           obj.iIdentifier = id;
       end
       
       function prevNode = get.hPrevNode(obj)
            prevNode = obj.hPrevNode;
       end
       
       function set.hPrevNode(obj, prevNode)
           
           obj.hPrevNode = prevNode;
       end
       
       function nextNode = get.hNextNode(obj)
            nextNode = obj.hNextNode;
       end
       
       function set.hNextNode(obj, nextNode)
           
           obj.hNextNode = nextNode;
       end
   end
    
end