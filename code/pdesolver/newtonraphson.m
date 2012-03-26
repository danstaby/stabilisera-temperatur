function [uGuess wGuess TGuess res n]  = newtonraphson(fu, fw, fT, uGuess, ...
					   wGuess,TGuess, uFreeNodes, ...
					   wFreeNodes, TFreeNodes, ...
					   uRange,wRange,TRange,ITERMAX, epsilon)

res = epsilon+1;
n = 0;
resvec = zeros(max(size(wGuess)+size(uGuess)+size(TGuess)),1);
x = zeros(max(size(resvec)),1);

while ((res > epsilon || n < 2) && n < ITERMAX)
  resvec(uRange,1) = LHSu(uGuess,wGuess,TGuess,uFreeNodes,wFreeNodes,TFreeNodes) ...
                   - fu(uFreeNodes);
  resvec(wRange,1) = LHSw(uGuess,wGuess,TGuess,uFreeNodes, wFreeNodes, ...
			  TFreeNodes) - fw(wFreeNodes);
  resvec(TRange,1) = LHST(uGuess,wGuess,TGuess,uFreeNodes,wFreeNodes,TFreeNodes)...
                   - fT(TFreeNodes);

  J = calcJacobian(uGuess, wGuess, TGuess, uFreeNodes, wFreeNodes, ...
		   TFreeNodes, uRange, wRange, TRange);

  dx = J\resvec;  
  x = x - dx;
  res = norm(resvec);
  n = n + 1;
  uGuess = x(uRange);
  wGuess = x(wRange);
  TGuess = x(TRange);

  fprintf(1, '%i\t\t%6.4g\t%6.4g\t%6.4g\n', n, norm(resvec(TRange)),...
          norm(resvec(uRange)), norm(resvec(wRange)))
end

