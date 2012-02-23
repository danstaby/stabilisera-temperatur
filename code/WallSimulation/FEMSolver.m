classdef FEMSolver < handle
    %Class that solves partial differential equations with the finite
    %element method
    %
    %TODO:
    %Add boundary conditions to GridNode
    %Add setting for weak formulation
    %Add boundary detection
    %Implement stiffness matrix and load vector generation
    %Implement solution presentation
    
    properties
        vGridSize
    end
    
    properties(GetAccess = private)
        hNodes
        hTriangles
    end
    
    methods
        function obj = FEMSolver(GridWidth, GridHeight)
           if(nargin < 1)
              GridWidth = 1;
              GridHeight = 1;
           end
            
           obj.vGridSize = [GridWidth GridHeight];
           obj.hNodes = LinkedList;
           obj.hTriangles = LinkedList;
        end
        
        
        function id = addTriangle(obj, nodes)
            hT = GridTriangle;
            
            
            for n=3:-1:1
               hN(n) = obj.hNodes.getElement(nodes(n)); 
            end
            
            hT.setNodes(nodes, hN);
            id = obj.hTriangles.addElement(hT);
        end
        
        function id = addNode(obj, pos, boundary)
            [id, hNode] = obj.addAndGetNode(pos, boundary);
        end
        
        function [id, hN] = addAndGetNode(obj, pos, boundary)
            id = obj.findNodeAtPosition(pos);
            if(id == 0)
                hN = GridNode;
                hN.vPosition = pos;

                if(nargin == 3)
                    hN.bBoundary = boundary;
                end

                id = obj.hNodes.addElement(hN);
            else
                hN = obj.hNodes.getElement(id);
            end
        end
        
        function ID = findNodeAtPosition(obj, pos)
            [hN, i] = obj.hNodes.getFirstElementAndID();
            ID = 0;
            while(i ~= 0)
                vPos = hN.vPosition;
                if(vPos(1) == pos(1) && vPos(2) == pos(2))
                    ID = i;
                    i = 0;
                else
                    [hN, i] = obj.hNodes.getNextElementAndID(i);
                end
                
            end
        end
        
        function splitTriangleMultipleTimes(obj, Identifier, times)
            %disp(['Splitting ' num2str(Identifier) ' multiple times'])
            if(times > 0)
               
               [idT, idN] = obj.splitTriangle(Identifier, 0);
               
               %disp(['Times: ' num2str(times)])
               for n=1:max(size(idT))
                  %idT(n)
                  obj.splitTriangleMultipleTimes(idT(n), times-1); 
               end
            end
            
        end
         
        function [idT, idN] = splitTriangle(obj, Identifier, visible)
            if(nargin < 3)
                visible = 1;
            end
            
           hOldT = obj.hTriangles.getElement(Identifier);
            
            
            nodes = hOldT.vNodes;
            pos = zeros(3,2);
            
            
            for n=3:-1:1
               hOldN(n) = obj.hNodes.getElement(nodes(n)); 
               pos(n,:) = hOldN(n).vPosition;
            end
           
            %Find mid points
            
            mp = zeros(3,2);
            mp(1,:) = 0.5*(pos(1,:)+pos(2,:));
            mp(2,:) = 0.5*(pos(1,:)+pos(3,:));
            mp(3,:) = 0.5*(pos(2,:)+pos(3,:));
            
            
            %Create new nodes and triangles
            hNewT(4) = GridTriangle;
            iNewT(4) = obj.hTriangles.addElement(hNewT(4));
            
            for n=3:-1:1
               [iNewN(n) hNewN(n)] = obj.addAndGetNode(mp(n,1:2));
               hNewT(n) = GridTriangle;
               iNewT(n) = obj.hTriangles.addElement(hNewT(n));
            end
            
            for n=1:4
               hNewT(n).bDraw = visible; 
            end
            
            %Set node points to triangles
            hNewT(1).setNodes([nodes(1) iNewN(1) iNewN(2)], ...
                              [hOldN(1) hNewN(1) hNewN(2)]);
                          
            hNewT(2).setNodes([nodes(2) iNewN(1) iNewN(3)], ...
                              [hOldN(2) hNewN(1) hNewN(3)]);
            
            hNewT(3).setNodes([nodes(3) iNewN(2) iNewN(3)], ...
                              [hOldN(3) hNewN(2) hNewN(3)]);
            
            hNewT(4).setNodes([iNewN(1) iNewN(2) iNewN(3)],...
                               [hNewN(1) hNewN(2) hNewN(3)]);
                           
            obj.removeTriangle(Identifier);
            
            idT = iNewT;
            idN = iNewN;
           
        end
        
        function Identifier = getTriangleIDFromPosition(obj, pos)
            Identifier = 0;
            [hT, ID] = obj.hTriangles.getFirstElementAndID();
            
            while(ID ~= 0)
                if(hT.isInside(pos) == 1)
                    Identifier = ID;
                    ID = 0;
                else
                    [hT, ID] = obj.hTriangles.getNextElementAndID(ID);
                    %disp(['Checking ID: ', num2str(ID)])
                end
            end
        end
        
        function pos = getNodePosition(obj, id)
            hN = obj.hNodes.getElement(id);
            
            if(hN ~= 0) 
                pos = hN.vPosition;
            else
                pos = 0;
            end
        end
        
        function pos = getTriangleCorners(obj, id)
            hT = obj.hTriangles.getElement(id);
            
            if(hT ~= 0)
                nodes = hT.hNodes;
                pos = zeros(3,2);
                for n=1:3
                    pos(n,:) = nodes(n).vPosition;
                end
            else
                pos = 0;
            end
        end
        
        function clearGrid(obj)
            obj.hNodes.clearList();
            obj.hTriangles.clearList();
            
        end
        
        function val = getNodeCount(obj)
            val = obj.hNodes.iSize;
        end
        
        function val = getTriangleCount(obj)
            val = obj.hTriangles.iSize;
        end
        
        function removeNode(obj, id)
            %Remove node and recalculate boundary
            
            %Find all triangles that contain the node
            %Remove those triangles
            %If the removed node lies on the boundary create new boundary
            %else fill the void with new triangles
            %
        end
    end
    
    methods(Access=private)
        
        function removeTriangle(obj, id)
            obj.hTriangles.removeElement(id);
        end
        
    end
    
end