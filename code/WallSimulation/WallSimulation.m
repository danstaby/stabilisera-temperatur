classdef WallSimulation < handle
    %Main worker for WallSimulation
    
    properties(GetAccess=private)
        %Calculation parameters
        
        %GUI handles
        hGUI
        hCmdRun
        hCmdLoadBoundary
        hTxtTriangleCount
        hTxtNodeCount
        hAxes
        
        %Object handles
        hFEM
        hNodeCircles
        hTriangleLines
        hSettings
        
        
        %Constants
        PlotHigh
        PlotLow
    end
    
    methods
        function obj = WallSimulation()
            %Constructor
            
            obj.PlotHigh = 1;
            obj.PlotLow = 0;
            
            obj.hFEM = FEMSolver;
            obj.hNodeCircles = LinkedList;
            obj.hTriangleLines = LinkedList;
            
            %Create graphical elements
            obj.hGUI = figure('Visible','off','Position',[360,500,650,385], ...
           'WindowStyle','modal', 'Name', 'Wall Simulation', ...
           'NumberTitle', 'off'); 

            obj.hCmdRun    = uicontrol('Style','pushbutton', ...
                         'String','Run Simulation','Position',[500,50,100,25]);

            obj.hCmdLoadBoundary    = uicontrol('Style','pushbutton', ...
                         'String','Load Boundary','Position',[500,80,100,25]);

            obj.hTxtTriangleCount = uicontrol('Style', 'text', ...
                                    'Position', [500,140,100,25], ...
                                     'BackgroundColor', get(obj.hGUI,'Color'));
            obj.hTxtNodeCount = uicontrol('Style', 'text', ...
                                    'Position', [500,170,100,25], ...
                                     'BackgroundColor', get(obj.hGUI,'Color'));

            obj.hAxes = axes('Units','Pixels','Position',[50,30,350,300],...
                              'XTick', [],'YTick', [], 'XLimMode' , 'manual',...
                              'YLimMode' , 'manual');
            axis(obj.hAxes, [obj.PlotLow obj.PlotHigh obj.PlotLow obj.PlotHigh])
            set(obj.hCmdLoadBoundary, 'Callback', @obj.promtLoadBoundary);
            set(obj.hAxes, 'ButtonDownFcn', @obj.divideTriangle);
            %obj.setupNewGrid(); %setup new grid
            obj.redrawGrid();
                          
            align([obj.hCmdRun],'Center','None');

            set(obj.hGUI, 'Visible', 'on') %Display form 

            
        end
        
        function promtLoadBoundary(obj, hObject, eventdata, var1, var2)
            
            strThisFile = mfilename('fullpath');
            n = strfind(strThisFile, '/');
            n = max(n);
            strThisFile = strThisFile(1:n);
            [FileName, PathName, FilterIndex] = ...
                uigetfile([strThisFile 'config/*.xls'], 'Load boundary file');
            
            if(FilterIndex ~= 0) %Check if the "Open" button was pressed
                
                obj.hSettings = Settings;
                
                obj.hSettings.ReadFile([PathName FileName]);
                
                obj.importSettings();
                
            end
            
        end
        
        function importSettings(obj)
            obj.hFEM.clearGrid();
            
            cNodes = obj.hSettings.cNodes;
            cBound = obj.hSettings.cBoundary;
            cTriangles = obj.hSettings.cTriangles;
            obj.importNodes(cBound, 1);
            obj.importNodes(cNodes, 0);
            
            obj.importTriangles(cTriangles);
            obj.updateGridText();
            
            
        end
        
        function importTriangles(obj, C)
            for n=1:max(size(C))
               row = C{n};
               id = obj.hFEM.addTriangle([row{1} row{2} row{3}]);
               obj.paintTriangle(id);
            end
        end
        
        function importNodes(obj, C, bound)
            
            for n=1:max(size(C))
                row = C{n};
               
                id = obj.hFEM.addNode([row{1} row{2}],bound);
                obj.paintNode(id);
            end
        end
        
        function redrawGrid(obj)
            cla(obj.hAxes)
            line([obj.PlotLow obj.PlotHigh], ...
                [obj.PlotHigh obj.PlotHigh], 'Parent' ,obj.hAxes, 'Color', 'k');
            line([obj.PlotHigh obj.PlotHigh], [obj.PlotHigh obj.PlotLow], ...
                'Parent', obj.hAxes, 'Color' ,'k');
                
        end
        
        function divideTriangle(obj, hObject, eventdata, var1, var2)
            pt = get(hObject, 'CurrentPoint');
            ID = obj.hFEM.getTriangleIDFromPosition(pt(1,1:2));
            
            if(ID ~= 0)
               [idT, idN] = obj.hFEM.splitTriangle(ID);
                obj.updateGridText();
                
                for n = 1:3
                   obj.paintNode(idN(n));
                   obj.paintTriangle(idT(n));
                end
                
                obj.paintTriangle(idT(4));
                %Draw new lines and new nodes
            end
        end
        
%         function setupNewGrid(obj, hObject, eventdata, var1, var2)
%             obj.redrawGrid();
%             
%             obj.hNodeCircles.clearList();
%             obj.hTriangleLines.clearList();
%             obj.hFEM.clearGrid();
%       
%             obj.createNode([0 0]);
%             obj.createNode([1 0]);
%             obj.createNode([0.5 0.5]);
%             obj.createNode([0 1]);
%             obj.createNode([1 1]);
%             
%             obj.createTriangle([1 2 3]);
%             obj.createTriangle([1 3 4]);
%             obj.createTriangle([2 3 5]);
%             obj.createTriangle([3 4 5]);
%             
%             obj.updateGridText();
%             
%         end
        
        function paintNode(obj, ID)
            h = obj.hNodeCircles.getElement(ID);
            if(h==0)
                pos = obj.hFEM.getNodePosition(ID);

                h = rectangle('Curvature', [1 1], 'Position' , ...
                        [pos(1)-0.01 pos(2)-0.01 0.02 0.02], 'Parent', obj.hAxes, ...
                        'LineWidth', 1, 'DisplayName', num2str(ID), ...
                        'FaceColor', 'k'); %Add click event here

                obj.hNodeCircles.addElement(h, ID);
            end
            
        end
        
        function paintTriangle(obj, ID)
            X = obj.hFEM.getTriangleCorners(ID);
            
            a = line([X(1,1) X(2,1)], [X(1,2) X(2,2)],'Color', 'k',...
                    'Parent', obj.hAxes);
            b = line([X(1,1) X(3,1)], [X(1,2) X(3,2)],'Color', 'k',...
                    'Parent', obj.hAxes);
            c = line([X(2,1) X(3,1)], [X(2,2) X(3,2)],'Color', 'k',...
                    'Parent', obj.hAxes);
            
            obj.hTriangleLines.addElement([a b c], ID);
            
        end
        
        function ID = createNode(obj, pos, boundary)
            if(nargin == 2)
                boundary = 0;
            end
            
            ID = obj.hFEM.addNode(pos, boundary);
            obj.paintNode(ID);
        end
        
        function ID = createTriangle(obj, nodes)
           ID = obj.hFEM.addTriangle(nodes);
           obj.paintTriangle(ID);
        end
        
        function updateGridText(obj)
            set(obj.hTxtNodeCount,'String', ...
                ['Node Count: ' num2str(obj.hFEM.getNodeCount)]);
            set(obj.hTxtTriangleCount, 'String', ...
                ['Triangle Count: ' num2str(obj.hFEM.getTriangleCount)]);
        end
    end
end