classdef Settings < handle
    
   properties
       cBoundary
       cNeumann
       cDirichlet
       cNodes
       cTriangles
       cWeakForm
       
   end
   
   properties(GetAccess=private)
       strFileName
   end
   
   methods
       function obj = Settings(~, strFile)
          if(nargin == 2)
              obj.strFileName = strFile;
          end
       end
       
       function ReadFile(obj, strFile)
          if(nargin == 2)
              obj.strFileName = strFile;
          end
          
          hFileLoader = FileLoader(obj.strFileName);
          hFileLoader.ReadFile();
          obj.cBoundary = hFileLoader.getField('#BOUND');
          obj.cNeumann = hFileLoader.getField('#BCNEU');
          obj.cDirichlet = hFileLoader.getField('#BCDIR');
          obj.cNodes = hFileLoader.getField('#NODES');
          obj.cTriangles = hFileLoader.getField('#ZONES');
          obj.cWeakForm = hFileLoader.getField('#WEAK');
          
          clear hFileLoader;
       end
      
   end
end