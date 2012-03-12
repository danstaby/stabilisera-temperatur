function M = createBoundaryMatrix(dim, borderCount, dbound, nbound)

M = [];
cCols = cell(1,borderCount);
maxSize = 0;
for n = 1:borderCount
  dC = dbound{n};
  nC = nbound{n};
  tmpQ = '0'*ones(1,dim.^2);
  tmpH = '0'*ones(1,dim);

  tmpH(1:(dim+1):dim^2) = '1';
  
  tmpQl = ones(1,dim.^2);
  tmpGl = ones(1,dim);
  tmpHl = ones(1,dim);
  tmpRl = ones(1,dim);
  
  tmpR = [];
  tmpG = [];
  for k = 1:dim
    %prep setting vectors
    dircount = 0;
    if(isnan(dC{k}) == 0)
      %insert into vector
      tmpRl(k) = max(size(dC{k}));
      tmpR = [tmpR char(dC{k})];
      dircount = dircount + 1;
    else
      tmpRl(k) = 1;
      tmpR = [tmpR '0'];
      tmpH((k-1)*dim+k) = '0';
    end

    if(isnan(nC{k}) == 0)
      %insert into vector
      tmpGl(k) = max(size(nC{k}));
      tmpG = [tmpG char(nC{k})];
      tmpQ(k) = '0';
    else
      tmpGl(k) = 1;
      tmpG = [tmpG '0'];
    end
  end
  
  if(dircount > 0)
    col = [dim, dircount, tmpQl, tmpGl, tmpHl, tmpRl, tmpQ, tmpG, tmpH, tmpR]';
  else
    col = [dim, 0, tmpQl, tmpGl, tmpQ, tmpG];
  end

  if(maxSize < max(size(col)))
    maxSize = max(size(col));
  end

  cCols{n} = col;
end

M = zeros(maxSize, borderCount);

for n = 1:borderCount
  col = cCols{n};
  s = max(size(col));
  M(1:s,n) = col;
end