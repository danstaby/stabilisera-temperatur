classdef FileLoader < handle
    
    properties
        strFileName
        
    end
    
    properties(GetAccess=private)
       cData 
       vHeaders
    end
    
    methods
        function obj = FileLoader(strName)
            if(nargin == 1)
               obj.strFileName = strName;
            end
        end
        
        function ReadFile(obj, strName)
            if(nargin == 2)
                obj.strFileName = strName;
            end
            
            
           [~,~,raw] = xlsread(obj.strFileName);
           
           obj.parseFile(raw);
           
           
        end
        
        function field = getField(obj, strHeader)
            ID = obj.getFieldID(strHeader);
            
            if(ID == 0)
                field = NaN;
                return;
            end
            
            field = obj.cData{ID};
            
        end
        
        function ID = getFieldID(obj, strHeader)
            ID = 0;
            
            for n = 1:max(size(obj.vHeaders))
               
               if(strcmp(obj.vHeaders{n},strHeader) == 1)
                   ID = n;
                   break;
               end
            end
        end
    end
    
    methods(Access=private)
        
        function parseFile(obj,raw)
           headCount = FileLoader.countHeaders(raw);
           obj.vHeaders = cell(headCount,1);
           
           obj.cData = cell(headCount,1);
           
           s = max(size(raw));
           headID = 0;
           tmpCell = {};
           rowID = 0;
           for n=1:s
               str = raw(n,1);
               if(iscellstr(raw(n,1)) == 1)
                   str = char(str);
                   
                   if(str(1) == '#')
                       
                       if(headID > 0)
                           %Put old header fields into cell
                           obj.cData{headID} = tmpCell;
                           tmpCell = {};
                           rowID = 0;
                       end

                       headID = headID + 1;
                       obj.vHeaders(headID) = cellstr(str);

                   else
                       %Put row into tmp cell

                       rowID = rowID + 1;

                       sr = size(raw(n,:));
                       LastNaN = sr(2)+1;
                       for m=sr(2):-1:1
                           if(isnan(raw{n,m})==1)
                               LastNaN = m;
                           end
                       end
                       LastNaN
                       tmpCell{rowID} = raw(n,1:(LastNaN-1))
                   end
               else
                  %Put row into tmp cell

                   rowID = rowID + 1;

                   sr = size(raw(n,:));
                   LastNaN = sr(2)+1;
                   for m=sr(2):-1:1
                        if(isnan(raw{n,m})==1)
                            LastNaN = m;
                        end
                   end
                   tmpCell{rowID} = raw(n,1:(LastNaN-1));
               end
           end
            
        end
        
    end
    
    methods(Static)
        function count = countHeaders(C)
            count = 0;
            s = size(C);
            for n = 1:s(1)
                str = C(n,1);
                if(iscellstr(C(n,1)) == 1)
                    str = char(str);

                    if(str(1) == '#')
                        count = count + 1; 
                    end
                end
            end
        end
    end
end